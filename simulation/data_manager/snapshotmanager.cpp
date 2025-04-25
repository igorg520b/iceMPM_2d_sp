#include "snapshotmanager.h"
#include "model.h"
#include "poisson_disk_sampling.h"

#include <spdlog/spdlog.h>
#include <H5Cpp.h>

#include <filesystem>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <utility>
#include <type_traits>

#include <fmt/format.h>
#include <fmt/std.h>

//#include <format>


#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

namespace fs = std::filesystem;


void icy::SnapshotManager::PrepareGrid(std::string fileNamePNG, std::string fileNameModelledArea)
{
    LOGR("icy::SnapshotManager::PrepareGrid {}", fileNamePNG);
    LOGR("modeled area file: {}", fileNameModelledArea);

    // (1) Load from HDF5: attributes, normals and path_indices
    H5::H5File file(fileNameModelledArea, H5F_ACC_RDONLY);
    int nPaths;

    H5::DataSet ds_path_indices = file.openDataSet("path_indices");
    ds_path_indices.openAttribute("nPaths").read(H5::PredType::NATIVE_INT, &nPaths);
    ds_path_indices.openAttribute("width").read(H5::PredType::NATIVE_INT, &model->prms.InitializationImageSizeX);
    ds_path_indices.openAttribute("height").read(H5::PredType::NATIVE_INT, &model->prms.InitializationImageSizeY);

    ds_path_indices.openAttribute("ModeledRegionOffsetX").read(H5::PredType::NATIVE_INT, &model->prms.ModeledRegionOffsetX);
    ds_path_indices.openAttribute("ModeledRegionOffsetY").read(H5::PredType::NATIVE_INT, &model->prms.ModeledRegionOffsetY);
    ds_path_indices.openAttribute("GridXTotal").read(H5::PredType::NATIVE_INT, &model->prms.GridXTotal);
    ds_path_indices.openAttribute("GridYTotal").read(H5::PredType::NATIVE_INT, &model->prms.GridYTotal);

    const int &width = model->prms.InitializationImageSizeX;
    const int &height = model->prms.InitializationImageSizeY;
    const int &ox = model->prms.ModeledRegionOffsetX;
    const int &oy = model->prms.ModeledRegionOffsetY;
    const int &gx = model->prms.GridXTotal;
    const int &gy = model->prms.GridYTotal;
    model->prms.cellsize = model->prms.DimensionHorizontal / (model->prms.InitializationImageSizeX-1);
    model->prms.cellsize_inv = 1.0/model->prms.cellsize;


    LOGR("initialization image: {} x {}", model->prms.InitializationImageSizeX, model->prms.InitializationImageSizeY);
    LOGR("grid size: {} x {}", model->prms.GridXTotal, model->prms.GridYTotal);
    LOGR("modeled area offset: {}, {}", model->prms.ModeledRegionOffsetX, model->prms.ModeledRegionOffsetY);

    const size_t sz = width*height;

    // temporary buffers for HDF5 data
    std::vector<Eigen::Vector2f> normals(sz);
    std::vector<int> path_indices(sz);

    ds_path_indices.read(path_indices.data(), H5::PredType::NATIVE_INT);
    file.openDataSet("normals").read(normals.data(), H5::PredType::NATIVE_FLOAT);
    file.close();



    // (2) allocate grid arrays host-side
    model->gpu.allocate_host_arrays_grid();



    // (3) Load PNG image (only for rendering)
    int channels, imgx, imgy;
    unsigned char *png_data = stbi_load(fileNamePNG.c_str(), &imgx, &imgy, &channels, 3); // expect 3 channels - RGB
    if(!png_data || channels != 3 || imgx != width || imgy != height)
    {
        LOGR("filename {} not loaded; channels {}", fileNamePNG, channels);
        LOGR("png expected size {} x {}; actual size {} x {}", width, height, imgx, imgy);
        throw std::runtime_error("png1 not loaded");
    }

    // function for obtaining index in png_data from the pixel's 2D index (i,j)
    auto idxInPng = [&](int i, int j) -> int { return 3*((height - j - 1)*width + i); };

    // save original colors for the whole image
    for(int i=0;i<width;i++)
        for(int j=0;j<height;j++)
            for(int k=0;k<3;k++)
                model->gpu.original_image_colors_rgb[(i+j*width)*3+k] = png_data[idxInPng(i, j)+k];

    stbi_image_free(png_data);



    // (4) Set grid normals and "status" byte
    auto transformPathIdx = [](const int &idx) -> uint8_t {
        if(idx < 1000) return (uint8_t)(idx + 1);
        else return (uint8_t)(idx-900);
    };


    // transfer / set grid land bit
    for(int i=0;i<gx;i++)
        for(int j=0;j<gy;j++)
        {
            size_t idx1 = (i+ox) + (j+oy)*width;
            size_t idx_host = j + i*gy; // store as row-major
            model->gpu.grid_status_buffer[idx_host] = transformPathIdx(path_indices[idx1]);
            model->gpu.grid_boundary_normals[idx_host] = (t_GridReal)normals[idx1].x();
            model->gpu.grid_boundary_normals[idx_host + gx*gy] = (t_GridReal)normals[idx1].y();
        }

}



void icy::SnapshotManager::PopulatePoints(std::string fileNameModelledAreaHDF5)
{
    LOGR("icy::SnapshotManager::PopulatePoints: {}", fileNameModelledAreaHDF5);

    // (0) params
    const int &width = model->prms.InitializationImageSizeX;
    const int &height = model->prms.InitializationImageSizeY;
    const int &ox = model->prms.ModeledRegionOffsetX;
    const int &oy = model->prms.ModeledRegionOffsetY;
    const int &gx = model->prms.GridXTotal;
    const int &gy = model->prms.GridYTotal;
    const size_t sz = width*height;
    const t_PointReal &h = model->prms.cellsize;

    std::vector<uint8_t> iceStatus(sz);     // 0 for open water, 1 for crushed, 2 for "intact"
    std::vector<float> iceThickness(sz);    // quantitive measure of "thickness"

    // (1) load iceStatus from HDF5
    H5::H5File file(fileNameModelledAreaHDF5, H5F_ACC_RDONLY);
    file.openDataSet("iceStatus").read(iceStatus.data(), H5::PredType::NATIVE_UINT8);
    file.openDataSet("iceThickness").read(iceThickness.data(), H5::PredType::NATIVE_FLOAT);
    file.close();


    // (2) generate points (or load from cache)
    std::vector<std::array<float, 2>> pt_buffer;   // initial buffer for reading/generating points
    generate_points(gx, gy, SimParams::MPM_points_per_cell, pt_buffer);
    const size_t initial_pt_count = pt_buffer.size();
    model->prms.ParticleVolume = h*h*gx*gy/ initial_pt_count;


    // convert unscaled point coordinates to (i,j) pair on the grid
    auto idxPt = [&](const std::array<float,2> &pt) -> std::pair<int,int> {
        const double scale = gx-1;
        return {(int)(pt[0]*scale + 0.5), (int)(pt[1]*scale + 0.5)};
    };


    // return true if the underlying cell is land or open water
    auto shouldRemove = [&](const std::array<float,2> &pt) -> bool {
        auto [i,j] = idxPt(pt);
        if(i<=1 || j<=1 || i>= (gx-2) || j>= (gy-2)) return true; // exclude points right at the boundary

        uint8_t status = iceStatus[(i+ox) + (j+oy)*width];
        return (status != 1 && status != 2);    // 1 - crushed, 2 - intact ice
    };


    std::erase_if(pt_buffer, shouldRemove);  // C++ 20
    model->prms.nPtsInitial = pt_buffer.size();


    // (3) allocate host-side array for points
    // allocate space for points on hssoa
    model->gpu.allocate_host_arrays_points();

    model->gpu.hssoa.size = model->prms.nPtsInitial; // set the actual number of points stored in HSSOA
    LOGR("icy::SnapshotManager::PopulatePoints: allocated {} points on HSSOA", model->prms.nPtsInitial);


    // (6) transfer points
    const double pointScale = (gx-1) * h;

    HostSideSOA &hssoa = model->gpu.hssoa;
    for(size_t k = 0; k<pt_buffer.size(); k++)
    {
        std::array<float,2> &pt = pt_buffer[k];

        auto [i,j] = idxPt(pt);

        // write into SOA
        SOAIterator it = hssoa.begin()+k;
        ProxyPoint &p = *it;
        p.setValue(SimParams::posx, pt[0]*pointScale);      // set point's x-position in SOA
        p.setValue(SimParams::posx+1, pt[1]*pointScale);    // set point's y-position in SOA

        const size_t idx_in_image = ((i+ox) + (j+oy)*width)*3;
        uint32_t r = model->gpu.original_image_colors_rgb[idx_in_image+0];
        uint32_t g = model->gpu.original_image_colors_rgb[idx_in_image+1];
        uint32_t b = model->gpu.original_image_colors_rgb[idx_in_image+2];
        uint32_t rgba = (r << 16) | (g << 8) | b;
        model->gpu.point_colors_rgb[k] = rgba;
        p.setValueInt(SimParams::integer_point_idx, k);

        // set ice type
        int idx = (i+ox) + (j+oy)*width;
        float thickness = iceThickness[idx];
        p.setValue(SimParams::idx_thickness, thickness);

        uint8_t status = iceStatus[idx]; // 1 - crushed; 2 - solid
        uint32_t val = 0;
        if(status == 1) val |= 0x10000;     // crushed
        p.setValueInt(SimParams::idx_utility_data, val);

        p.setValue(SimParams::idx_Jp_inv, 1.f);
        p.setValue(SimParams::idx_Qp, 1.f);
        for(int idx=0; idx<SimParams::dim; idx++)
            p.setValue(SimParams::Fe00+idx*2+idx, 1.f);
    }
    LOGV("Points transferred to HSSOA");

    // convert points to cell-based local coordinates
    hssoa.convertToIntegerCellFormat(h);

    // we no longer need the original color in the sattelite photo
    FillModelledAreaWithBlueColor();

    if(model->prms.SaveSnapshots)
    {
        SavePointColors();
        SaveSnapshot(0, 0);
    }
}



void icy::SnapshotManager::FillModelledAreaWithBlueColor()
{
    const int &width = model->prms.InitializationImageSizeX;
    const int &height = model->prms.InitializationImageSizeY;
    const int &ox = model->prms.ModeledRegionOffsetX;
    const int &oy = model->prms.ModeledRegionOffsetY;
    const int &gx = model->prms.GridXTotal;
    const int &gy = model->prms.GridYTotal;

    for(int i=0;i<gx;i++)
        for(int j=0;j<gy;j++)
        {
            uint8_t status = model->gpu.grid_status_buffer[j + i*gy];
            if(status == 100)
            {
                for(int k=0;k<3;k++)
                    model->gpu.original_image_colors_rgb[((i+ox)+(j+oy)*width)*3+k] = waterColor[k];
            }
        }
}

void icy::SnapshotManager::SplitIntoPartitionsAndTransferToDevice()
{
    // particle volume and mass
    model->prms.ComputeHelperVariables();

    // allocate GPU partitions
    model->gpu.initialize();
    model->gpu.split_hssoa_into_partitions();
    model->gpu.allocate_arrays();
    model->gpu.transfer_to_device();
}








void icy::SnapshotManager::ReadPointsFromSnapshot(std::string fileNameSnapshotHDF5)
{
    LOGR("ReadPointsFromSnapshot {}", fileNameSnapshotHDF5);
    // instead of generating the points, read them from the provided HDF5 (snapshot) file
    H5::H5File file(fileNameSnapshotHDF5, H5F_ACC_RDONLY);
    H5::DataSet ds = file.openDataSet("pts_data");

    ds.openAttribute("SimulationStep").read(H5::PredType::NATIVE_INT, &model->prms.SimulationStep);
    ds.openAttribute("SimulationTime").read(H5::PredType::NATIVE_DOUBLE, &model->prms.SimulationTime);
    ds.openAttribute("HSSOA_size").read(H5::PredType::NATIVE_UINT, &model->gpu.hssoa.size);

    int nPtsArrays;
    ds.openAttribute("nPtsArrays").read(H5::PredType::NATIVE_INT, &nPtsArrays);

    H5::DataSpace dsp = ds.getSpace();
    hsize_t dims[2];
    dsp.getSimpleExtentDims(dims, nullptr);

    if(dims[0] != SimParams::nPtsArrays || dims[1] != model->gpu.hssoa.size || nPtsArrays != SimParams::nPtsArrays)
        throw std::runtime_error("icy::SnapshotManager::ReadSnapshot array size mismatch");


    model->prms.nPtsInitial = model->gpu.hssoa.size;
    model->gpu.allocate_host_arrays_points();
    model->gpu.hssoa.size = model->prms.nPtsInitial;

    // Define hyperslab
    hsize_t dims_mem[2] = {SimParams::nPtsArrays, model->gpu.hssoa.capacity};
    H5::DataSpace memspace(2, dims_mem);
    hsize_t offset[2] = {0, 0};      // Start reading at the origin
    hsize_t count[2] = {dims[0], dims[1]}; // Size of the data in the file
    memspace.selectHyperslab(H5S_SELECT_SET, count, offset);

    H5::DataType dtype;
    if constexpr(std::is_same_v<t_PointReal, float>) dtype = H5::PredType::NATIVE_FLOAT;
    else dtype = H5::PredType::NATIVE_DOUBLE;

    ds.read(model->gpu.hssoa.host_buffer, dtype, memspace, dsp);


//    model->gpu.hssoa.RemoveDisabledAndSort(model->prms.GridYTotal);

    FillModelledAreaWithBlueColor();
    ReadPointColors();
    LOGR("ReadPointsFromSnapshot; hssoa capacity {}; size {}", model->gpu.hssoa.capacity, model->gpu.hssoa.size);
}


void icy::SnapshotManager::SaveSnapshot(int SimulationStep, double SimulationTime)
{
    // ensure that the output directory exists
    fs::path outputDir = "output";
    fs::path snapshotsDir = "snapshots";
    fs::path targetPath = outputDir / SimulationTitle / snapshotsDir;
    fs::create_directories(targetPath);

    // save current state
    const int frame = SimulationStep / model->prms.UpdateEveryNthStep;
//    std::string baseName = std::format("s{:05d}.h5", frame);
    std::string baseName = fmt::format(fmt::runtime("s{:05d}.h5"), frame);

    fs::path fullPath = targetPath / baseName;


    H5::H5File file(fullPath.string(), H5F_ACC_TRUNC);

    // points
    hsize_t dims_points[2] = {SimParams::nPtsArrays, model->gpu.hssoa.size};
    H5::DataSpace dataspace_points(2, dims_points);

    // We are zipping the data, so using chunks
    H5::DSetCreatPropList proplist;
    hsize_t chunk_size = (hsize_t)std::min((unsigned)256*1024, model->gpu.hssoa.size);
    hsize_t chunk_dims[2] = {SimParams::nPtsArrays, chunk_size};
    proplist.setChunk(2, chunk_dims);
    proplist.setDeflate(5);

    H5::DataType dtype;
    if constexpr(std::is_same_v<t_PointReal, float>) dtype = H5::PredType::NATIVE_FLOAT;
    else dtype = H5::PredType::NATIVE_DOUBLE;

    // Define the hyperslab in the memory space
    hsize_t mem_dims[2] = {SimParams::nPtsArrays, model->gpu.hssoa.capacity};
    H5::DataSpace memspace(2, mem_dims);
    hsize_t mem_offset[2] = {0, 0};  // Start at the beginning of the memory array
    hsize_t mem_count[2] = {SimParams::nPtsArrays, model->gpu.hssoa.size};  // Number of elements to select
    memspace.selectHyperslab(H5S_SELECT_SET, mem_count, mem_offset);

    // Write the data to the dataset
    H5::DataSet dataset_pts = file.createDataSet("pts_data", dtype, dataspace_points, proplist);
    dataset_pts.write(model->gpu.hssoa.host_buffer, dtype, memspace, dataspace_points);

    H5::DataSpace att_dspace(H5S_SCALAR);
    dataset_pts.createAttribute("SimulationStep", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &SimulationStep);
    dataset_pts.createAttribute("SimulationTime", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &SimulationTime);
    dataset_pts.createAttribute("HSSOA_size", H5::PredType::NATIVE_UINT, att_dspace).write(H5::PredType::NATIVE_UINT, &model->gpu.hssoa.size);
    int nPtsArrays = SimParams::nPtsArrays;
    dataset_pts.createAttribute("nPtsArrays", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &nPtsArrays);
}



void icy::SnapshotManager::SavePointColors()
{
    // ensure that the output directory exists
    fs::path outputDir = "output";
    fs::path snapshotsDir = "snapshots";
    fs::path targetPath = outputDir / SimulationTitle / snapshotsDir;
    fs::create_directories(targetPath);

    // save current state
    fs::path fullPath = targetPath / "point_colors.h5";

    H5::H5File file(fullPath.string(), H5F_ACC_TRUNC);

    hsize_t dims_points[1] = {model->gpu.point_colors_rgb.size()};
    H5::DataSpace dataspace_points(1, dims_points);

    // We are zipping the data, so using chunks
    H5::DSetCreatPropList proplist;
    hsize_t chunk_size = (hsize_t)std::min((hsize_t)256*1024, dims_points[0]);
    hsize_t chunk_dims[1] = {chunk_size};
    proplist.setChunk(1, chunk_dims);
    proplist.setDeflate(5);

    H5::DataSet dataset_pts = file.createDataSet("pts_colors", H5::PredType::NATIVE_UINT32, dataspace_points, proplist);
    dataset_pts.write(model->gpu.point_colors_rgb.data(), H5::PredType::NATIVE_UINT32);
}


void icy::SnapshotManager::ReadPointColors()
{
    // ensure that the output directory exists
    fs::path outputDir = "output";
    fs::path snapshotsDir = "snapshots";
    fs::path targetPath = outputDir / SimulationTitle / snapshotsDir;
    fs::create_directories(targetPath);

    // save current state
    fs::path fullPath = targetPath / "point_colors.h5";

    H5::H5File file(fullPath.string(), H5F_ACC_RDONLY);
    file.openDataSet("pts_colors").read(model->gpu.point_colors_rgb.data(), H5::PredType::NATIVE_UINT32);
}








std::string icy::SnapshotManager::prepare_file_name(int gx, int gy)
{
//    return std::format("{}/point_cache_{:05d}_{:05d}.h5", pts_cache_path, gx, gy);
    return fmt::format("{}/point_cache_{:05d}_{:05d}.h5", pts_cache_path, gx, gy);
}





// Haversine formula implementation
double icy::SnapshotManager::haversineDistance(double lat, double lon1, double lon2) {
    // Earth's radius in meters
    constexpr double R = SimParams::Earth_Radius;

    // Convert coordinates to radians
    double latRad = degreesToRadians(lat);
    double lon1Rad = degreesToRadians(lon1);
    double lon2Rad = degreesToRadians(lon2);

    // Difference in longitude
    double deltaLon = lon2Rad - lon1Rad;

    // Haversine formula
    double a = std::pow(std::sin(deltaLon / 2), 2);
    double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));

    // Distance along a parallel (latitude constant)
    double distance = R * c * std::cos(latRad);

    return distance;
}


/*
void icy::SnapshotManager::LoadWindData(std::string fileName)
{
    LOGR("icy::SnapshotManager::LoadWindData - {}", fileName);
    model->prms.use_GFS_wind = true;
    model->wind_interpolator.LoadWindData(fileName,
                                          model->prms.LatMin,model->prms.LatMax,model->prms.LonMin,model->prms.LonMax,
                                          model->prms.SimulationStartUnixTime);
    model->prms.gridLatMin = model->wind_interpolator.gridLatMin;
    model->prms.gridLonMin = model->wind_interpolator.gridLonMin;
}
*/



// =========================  POINTS GENERATION

void icy::SnapshotManager::generate_points(int gx, int gy, float points_per_cell,
                                                  std::vector<std::array<float, 2>> &buffer)
{
    bool cache_result = attempt_to_fill_from_cache(gx, gy, buffer);
    if(!cache_result) generate_and_save(gx, gy, SimParams::MPM_points_per_cell, buffer);
}


void icy::SnapshotManager::generate_and_save(int gx, int gy, float points_per_cell, std::vector<std::array<float, 2>> &buffer)
{
    const float dy = 1.0f*gy/gx;

    const std::array<float, 2>kXMin {0, 0};
    const std::array<float, 2>kXMax {1, dy};
    constexpr float magic_constant = 0.6;
    const float kRadius = sqrt(magic_constant/(points_per_cell*gx*gx));

    LOGR("generate_and_save: attempting to generate {} points in {}x{} grid", (points_per_cell*gx*gy), gx, gy);
    buffer = thinks::PoissonDiskSampling(kRadius, kXMin, kXMax);

    const float result_ppc = (float)buffer.size()/(gx*gy);
    LOGR("grid: {}x{}; generated pts {:>8}; density {:>7.4}",gx, gy, buffer.size(), result_ppc);

    const float scale = sqrt(result_ppc/(points_per_cell*1.0005));
    if(scale<1.0)
    {
        LOGR("requested ppc {}; generated ppc {}", points_per_cell, result_ppc);
        throw std::runtime_error("point generation error");
    }

    LOGR("requested ppc {}; generated ppc {}; scale {}%", points_per_cell, result_ppc, 100*(scale-1.f));

    for(auto &pt : buffer)
    {
        pt[0] *= scale;
        pt[1] *= scale;
    }

    auto result_it = std::remove_if(buffer.begin(),buffer.end(), [&](std::array<float,2> &pt)
                                    {return (pt[0]>1.f || pt[1]>dy ||pt[0] < 0 || pt[1] < 0);
                                    });
    buffer.erase(result_it, buffer.end());
    size_t final_count = buffer.size();
    float final_ppc = (float)final_count/(gx*gy);
    LOGR("grid: {}x{}; updated pts {:>8}; updated_density {:>7.4}",gx, gy, buffer.size(), final_ppc);

    // sort
    double hinv = (float)(gx-1);
    auto computeIdx = [&](std::array<float,2> &pt) -> int {
        float &x = pt[0];
        float &y = pt[1];

        int x_idx = (uint32_t)(x*hinv + 0.5);
        int y_idx = (uint32_t)(y*hinv + 0.5);
        int result_idx = x_idx*gy + y_idx;
        return result_idx;
    };

    LOGV("sorting started");
    std::sort(buffer.begin(),buffer.end(),
              [&](std::array<float,2> &pt1, std::array<float,2> &pt2)
              {return computeIdx(pt1) < computeIdx(pt2);}
              );

    LOGR("sorting finished; first point: {}, {}", buffer.front()[0], buffer.front()[1]);

    std::vector<int> grid_count(gx*gy,0);
    for(auto &pt : buffer) grid_count[computeIdx(pt)]++;

    auto [it_min, it_max] = std::minmax_element(grid_count.begin(),grid_count.end());

    LOGR("grid fill min/max {} to {} ", *it_min, *it_max);

    std::vector<int> histogram(std::max((int)points_per_cell*3, *it_max), 0);
    for(int &count : grid_count) histogram[count]++;
    for(auto val : histogram) std::cout << val << ' ';
    std::cout << std::endl;

    // write file
    std::filesystem::path directory_path(pts_cache_path);
    if (!std::filesystem::exists(directory_path)) std::filesystem::create_directories(directory_path);
    std::string fileName = prepare_file_name(gx, gy);
    LOGR("saving file {}", fileName);
    H5::H5File file(fileName, H5F_ACC_TRUNC);

    const size_t nPts = buffer.size();
    hsize_t dims_pts[2] = {nPts,2};
    H5::DataSpace dataspace_pts(2, dims_pts);

    H5::DataSet dataset = file.createDataSet("coords", H5::PredType::NATIVE_FLOAT, dataspace_pts);
    dataset.write(buffer.data(), H5::PredType::NATIVE_FLOAT);

    H5::DataSpace att_dspace(H5S_SCALAR);
    dataset.createAttribute("gx", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &gx);
    dataset.createAttribute("gy", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &gy);
    file.close();
}


bool icy::SnapshotManager::attempt_to_fill_from_cache(int gx, int gy, std::vector<std::array<float, 2>> &buffer)
{
    LOGV("attempting to load from cache");
    std::string fileName = prepare_file_name(gx, gy);
    std::filesystem::path file_path(fileName);
    if (!std::filesystem::exists(file_path))
    {
        LOGV("cached file does not exist");
        return false;
    }

    H5::H5File file(fileName, H5F_ACC_RDONLY);
    H5::DataSet dataset = file.openDataSet("coords");

    hsize_t dims[2];
    dataset.getSpace().getSimpleExtentDims(dims,NULL);
    if(dims[1] != 2) {
        LOGR("something is wrong with cache file - dims[1] is {}",dims[1]);
        return false;
    }

    int saved_gx, saved_gy;
    dataset.openAttribute("gx").read(H5::PredType::NATIVE_INT, &saved_gx);
    dataset.openAttribute("gy").read(H5::PredType::NATIVE_INT, &saved_gy);
    if(saved_gx != gx || saved_gy != gy)
    {
        LOGR("cache error - expected grid {}x{}; cached grid {}x{}",gx, gy, saved_gx, saved_gy);
        return false;
    }

    buffer.resize(dims[0]);
    LOGR("reading {} points from file", dims[0]);
    dataset.read(buffer.data(), H5::PredType::NATIVE_FLOAT);
    LOGR("attempt_to_fill_from_cache: finished reading from file; {} points", buffer.size());
    LOGR("attempt_to_fill_from_cache: first point: {}, {}", buffer.front()[0], buffer.front()[1]);

    return true;
}






// =============================  READ AND WRITE SNAPSHOTS

void icy::SnapshotManager::PrepareFrameArrays()
{
    const int gridSize = model->prms.GridXTotal*model->prms.GridYTotal;
    count.assign(gridSize,0); // in case we render a 3-channel image
    vis_vx.assign(gridSize,0);
    vis_vy.assign(gridSize,0);
    vis_r.assign(gridSize,0);
    vis_g.assign(gridSize,0);
    vis_b.assign(gridSize,0);
    vis_Jpinv.assign(gridSize,0);
    vis_P.assign(gridSize,0);
    vis_Q.assign(gridSize,0);
    vis_alpha.assign(gridSize,0);
    rgb.assign(gridSize*3,0);

    const int nPts = model->gpu.hssoa.size;

    for(int i=0;i<nPts;i++)
    {
        SOAIterator s = model->gpu.hssoa.begin()+i;
        //PointVector2r pos = s->getPos(model->prms.cellsize);
        int cellIdx = s->getCellIndex(model->prms.GridYTotal);

        count[cellIdx]++;
        vis_vx[cellIdx] += static_cast<float>(s->getValue(SimParams::velx+0));
        vis_vy[cellIdx] += static_cast<float>(s->getValue(SimParams::velx+1));

        int pt_idx = s->getValueInt(SimParams::integer_point_idx);
        uint32_t rgb = model->gpu.point_colors_rgb[pt_idx];
        uint8_t r = (rgb >> 16) & 0xff;
        uint8_t g = (rgb >> 8) & 0xff;
        uint8_t b = rgb & 0xff;

        vis_r[cellIdx] += r/255.;
        vis_g[cellIdx] += g/255.;
        vis_b[cellIdx] += b/255.;

        vis_Jpinv[cellIdx] += static_cast<float>(s->getValue(SimParams::idx_Jp_inv));
        vis_P[cellIdx] += static_cast<float>(s->getValue(SimParams::idx_P));
        vis_Q[cellIdx] += static_cast<float>(s->getValue(SimParams::idx_Q));
    }

#pragma omp parallel for
    for(int i=0;i<gridSize;i++)
    {
        int n = count[i];
        if(n>0)
        {
            float coeff = 1./n;
            vis_vx[i] *= coeff;
            vis_vy[i] *= coeff;
            vis_r[i] *= coeff;
            vis_g[i] *= coeff;
            vis_b[i] *= coeff;
            vis_Jpinv[i] *= coeff;
            vis_P[i] *= coeff;
            vis_Q[i] *= coeff;

            // Scale the float values to the range [0, 255]
            uint8_t red = static_cast<uint8_t>(vis_r[i] * 255.0f);
            uint8_t green = static_cast<uint8_t>(vis_g[i] * 255.0f);
            uint8_t blue = static_cast<uint8_t>(vis_b[i] * 255.0f);
            rgb[i*3+0] = red;
            rgb[i*3+1] = green;
            rgb[i*3+2] = blue;
        }
        else
        {
            vis_vx[i] = vis_vy[i] = vis_Jpinv[i] = vis_P[i] = vis_Q[i] = 0;
            rgb[i*3+0] = rgb[i*3+1] = rgb[i*3+2] = 0;
        }
    }

}



void icy::SnapshotManager::SaveFrame(int SimulationStep, double SimulationTime)
{
    LOGR("SnapshotManager::SaveFrame: step {}, time {}", SimulationStep, SimulationTime);
    PrepareFrameArrays();

    // save as HDF5
    fs::path outputDir = "output";
    fs::path framesDir = "frames";
    fs::path targetPath = outputDir / SimulationTitle / framesDir;
    fs::create_directories(targetPath);

    // save current state
    const int frame = SimulationStep / model->prms.UpdateEveryNthStep;
//    std::string baseName = std::format("f{:05d}.h5", frame);
    std::string baseName = fmt::format(fmt::runtime("f{:05d}.h5"), frame);

    fs::path fullPath = targetPath / baseName;
    H5::H5File file(fullPath.string(), H5F_ACC_TRUNC);

    hsize_t dims_grid[2] = {(hsize_t)model->prms.GridXTotal, (hsize_t)model->prms.GridYTotal};
    H5::DataSpace dataspace_grid(2, dims_grid);

    H5::DSetCreatPropList proplist;
    hsize_t chunk_dims[2] = {128, 128};
    proplist.setChunk(2, chunk_dims);
    proplist.setDeflate(6);

    H5::DataSet ds_vis_vx = file.createDataSet("vis_vx", H5::PredType::NATIVE_FLOAT, dataspace_grid, proplist);
    ds_vis_vx.write(vis_vx.data(), H5::PredType::NATIVE_FLOAT);

    H5::DataSet ds_vis_vy = file.createDataSet("vis_vy", H5::PredType::NATIVE_FLOAT, dataspace_grid, proplist);
    ds_vis_vy.write(vis_vy.data(), H5::PredType::NATIVE_FLOAT);

    H5::DataSet ds_vis_JpInv = file.createDataSet("vis_Jpinv", H5::PredType::NATIVE_FLOAT, dataspace_grid, proplist);
    ds_vis_JpInv.write(vis_Jpinv.data(), H5::PredType::NATIVE_FLOAT);

    H5::DataSet ds_vis_P = file.createDataSet("vis_P", H5::PredType::NATIVE_FLOAT, dataspace_grid, proplist);
    ds_vis_P.write(vis_P.data(), H5::PredType::NATIVE_FLOAT);

    H5::DataSet ds_vis_Q = file.createDataSet("vis_Q", H5::PredType::NATIVE_FLOAT, dataspace_grid, proplist);
    ds_vis_Q.write(vis_Q.data(), H5::PredType::NATIVE_FLOAT);

    H5::DataSet ds_count = file.createDataSet("count", H5::PredType::NATIVE_UINT8, dataspace_grid, proplist);
    ds_count.write(count.data(), H5::PredType::NATIVE_UINT8);

    hsize_t dims_rgb[3] = {(hsize_t)model->prms.GridXTotal, (hsize_t)model->prms.GridYTotal, 3};
    H5::DataSpace dataspace_rgb(3, dims_rgb);

    H5::DSetCreatPropList proplist2;
    hsize_t chunk_dims2[3] = {128, 128, 3};
    proplist2.setChunk(3, chunk_dims2);
    proplist2.setDeflate(6);

    H5::DataSet ds_rgb = file.createDataSet("rgb", H5::PredType::NATIVE_UINT8, dataspace_rgb, proplist2);
    ds_rgb.write(rgb.data(), H5::PredType::NATIVE_UINT8);

    H5::DataSpace att_dspace(H5S_SCALAR);
    ds_count.createAttribute("SimulationStep", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &SimulationStep);
    ds_count.createAttribute("SimulationTime", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &SimulationTime);
}




