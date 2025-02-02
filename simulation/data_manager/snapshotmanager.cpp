#include "snapshotmanager.h"
#include "model.h"
#include "poisson_disk_sampling.h"


#include <spdlog/spdlog.h>
#include <fmt/format.h>
#include <H5Cpp.h>

#include <filesystem>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <utility>
#include <type_traits>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


void icy::SnapshotManager::PreparePointsAndSetupGrid(std::string fileName, std::string fileNameModelledArea)
{
    spdlog::info("icy::SnapshotManager::PreparePointsAndSetupGrid {}",fileName);
    spdlog::info("modeled area file: {}", fileNameModelledArea);
//    unsigned char* png_data, *png_data_modelled_region;

    // load satellite photo
    int channels, imgx, imgy;
    unsigned char *png_data = stbi_load(fileName.c_str(), &imgx, &imgy, &channels, 4); // request 4 channels - RGBA
    if(!png_data || channels != 4)
    {
        spdlog::error("filename {}", fileName);
        spdlog::error("png_data {}", (void*)png_data);
        spdlog::error("channels {}", channels);
        spdlog::error("png1 not loaded");
        throw std::runtime_error("png1 not loaded");
    }
    model->prms.InitializationImageSizeX = imgx;
    model->prms.InitializationImageSizeY = imgy;

    // load highlighted modelling area
    unsigned char *png_data_modelled_region = stbi_load(fileNameModelledArea.c_str(), &imgx, &imgy, &channels, 4); // request 4 channels - RGBA
    if(imgx != model->prms.InitializationImageSizeX || imgy != model->prms.InitializationImageSizeY
        || !png_data_modelled_region || channels != 4)
    {
        spdlog::error("png2 not loaded: {}", fileNameModelledArea);
        throw std::runtime_error("png2 not loaded");
    }

    // function for obtaining index in png_data from the pixel's 2D index (i,j)
    auto idxInPng = [&](int i, int j) -> int { return 4*((imgy - j - 1)*imgx + i); };

    // function that determines whether (i,j) point at land in png_data_modelled_region
    auto is_modeled_area = [&](int i, int j) -> bool {
        unsigned char _a = png_data_modelled_region[idxInPng(i,j)+3];
        unsigned char _b = png_data_modelled_region[idxInPng(i,j)+2];
        return ((_a > 128) && (_b > 250));
    };

    // determine the extent and offset of the highlighted (modeled) region
    int xmin = imgx;
    int xmax = -1;
    int ymin = imgy;
    int ymax = -1;

    for(int i=0;i<imgx;i++)
        for(int j=0;j<imgy;j++)
        {
            if(is_modeled_area(i,j))
            {
                xmin = std::min(xmin, i);
                xmax = std::max(xmax, i);
                ymin = std::min(ymin, j);
                ymax = std::max(ymax, j);
            }
        }

    xmin = std::max(0, xmin-2);
    xmax = std::min(imgx-1, xmax+2);
    ymin = std::max(0, ymin-2);
    ymax = std::min(imgy-1, ymax+2);

    const int ox = model->prms.ModeledRegionOffsetX = xmin;
    const int oy = model->prms.ModeledRegionOffsetY = ymin;
    const int gx = model->prms.GridXTotal = xmax-xmin+1;
    const int gy = model->prms.GridYTotal = ymax-ymin+1;
    if(model->prms.GridXTotal < 0 || model->prms.GridYTotal < 0)
    {
        spdlog::error("modeled area not found");
        throw std::runtime_error("modeled area not found");
    }

    spdlog::info("initialization image: {} x {}", model->prms.InitializationImageSizeX, model->prms.InitializationImageSizeY);
    spdlog::info("grid size: {} x {}", model->prms.GridXTotal, model->prms.GridYTotal);
    spdlog::info("modeled area offset: {}, {}", model->prms.ModeledRegionOffsetX, model->prms.ModeledRegionOffsetY);
    spdlog::info("modeled area extent: [{}, {}] - [{}, {}]", xmin, ymin, xmax, ymax);

    // either load or generate points
    std::vector<std::array<float, 2>> buffer;   // initial buffer for reading/generating points
    generate_points(gx, gy, SimParams::MPM_points_per_cell, buffer);

    // convert unscaled point coordinates to (i,j) pair on the grid
    auto idxPt = [&](std::array<float,2> &pt) -> std::pair<int,int> {
        return {(int)((double)pt[0]*(gx-1) + 0.5), (int)((double)pt[1]*(gx-1) + 0.5)};
    };

    // now that we have png_data, classify points and count those remaining (not land nor open water)
    auto keepPoint = [&](std::array<float,2> &pt) -> bool {
        auto [i,j] = idxPt(pt);
        if(i<=1 || j<=1 || i>= (gx-2) || j>= (gy-2)) return false; // exclude points right at the boundary
        if(!is_modeled_area(i+ox, j+oy)) return false; // outside modelled region

        int idx_png = idxInPng(i+ox,j+oy);
        unsigned char r = png_data[idx_png + 0];
        unsigned char g = png_data[idx_png + 1];
        unsigned char b = png_data[idx_png + 2];
        Eigen::Vector3f rgb((float)r/255.,(float)g/255.,(float)b/255.);

        auto [category, interpValue] = categorizeColor(rgb);
        if(category <= 0) return false; // open water
        else return true; // some sort of ice (keep it)
    };

    model->prms.nPtsInitial = std::count_if(buffer.begin(), buffer.end(), keepPoint);

    // allocate host-side arrays for points and grid
    model->gpu.allocate_host_arrays();
    model->gpu.hssoa.size = model->prms.nPtsInitial; // set the actual number of points stored in HSSOA
    spdlog::info("PreparePointsAndSetupGrid: nPoints {}", model->prms.nPtsInitial);


    // compute the X-dimension from lat/lon using haversineDistance()
    const double XScale = model->prms.DimensionHorizontal;
    model->prms.cellsize = model->prms.DimensionHorizontal / (model->prms.InitializationImageSizeX-1);
    model->prms.cellsize_inv = 1.0/model->prms.cellsize;
    const t_PointReal &h = model->prms.cellsize;
    model->prms.ParticleVolume = h*h*gx*gy/ buffer.size();
    const double pointScale = (gx-1) * h;

    spdlog::info("assumed horizontal scale: {}; cellsize {}", XScale, model->prms.cellsize);


//    if(model->prms.DimensionHorizontal == 0)
//    {
//        XScale = haversineDistance((model->prms.LatMin + model->prms.LatMax)*0.5, model->prms.LonMin, model->prms.LonMax);
//    }
//    else
//    {
 //       XScale = model->prms.DimensionHorizontal;
 //   }

    // compute cellsize

    // sequentially transfer points from "buffer" to hssoa (only those counted), scale coordinates
    uint32_t count = 0;
    HostSideSOA &hssoa = model->gpu.hssoa;
    for(std::array<float,2> &pt : buffer)
    {
        if(!keepPoint(pt)) continue;
        if(count >= model->prms.nPtsInitial) throw std::runtime_error("when transferring from buffer to SOA - index error");

        auto [i,j] = idxPt(pt);
        int idx_png = idxInPng(i+ox, j+oy);

        unsigned char r = png_data[idx_png + 0];
        unsigned char g = png_data[idx_png + 1];
        unsigned char b = png_data[idx_png + 2];
        Eigen::Vector3f rgb((float)r/255.,(float)g/255.,(float)b/255.);
        auto [category, interpValue] = categorizeColor(rgb);
        if(category <= 0) continue; // either land or water (should not happen)

        // write into SOA
        SOAIterator it = hssoa.begin()+count;
        ProxyPoint &p = *it;
        p.setValue(SimParams::posx, pt[0]*pointScale);      // set point's x-position in SOA
        p.setValue(SimParams::posx+1, pt[1]*pointScale);    // set point's y-position in SOA

        uint32_t rgba = (static_cast<uint32_t>(r) << 16) |  (static_cast<uint32_t>(g) << 8) |  static_cast<uint32_t>(b);
        model->gpu.point_colors_rgb[count] = rgba;
        p.setValueInt(SimParams::integer_point_idx, count);

        uint32_t val = 0;
        if(category == 1)
        {
            val |= 0x10000;     // crushed
            p.setValueInt(SimParams::idx_utility_data, val);
            p.setValue(SimParams::idx_initial_strength, 0.3);
        }
        else
        {
            p.setValue(SimParams::idx_initial_strength, interpValue*0.5+0.5);
        }

        p.setValue(SimParams::idx_Jp_inv, 1.f);
        p.setValue(SimParams::idx_Qp, 1.f);
        for(int idx=0; idx<SimParams::dim; idx++)
            p.setValue(SimParams::Fe00+idx*2+idx, 1.f);

        count++;
    }

    // convert points to cell-based local coordinates
    hssoa.convertToIntegerCellFormat(h);


    // transfer / set grid land bit
    for(int i=0;i<gx;i++)
        for(int j=0;j<gy;j++)
        {
            bool is_modeled = is_modeled_area(i+ox,j+oy);
            int host_grid_idx = j + i*gy;
            model->gpu.grid_status_buffer[host_grid_idx] = is_modeled ? 1 : 0;
        }

    // save original colors for the whole image
    for(int i=0;i<imgx;i++)
        for(int j=0;j<imgy;j++)
        {
            for(int k=0;k<3;k++)
            {
                model->gpu.original_image_colors_rgb[(i+j*imgx)*3+k] = png_data[idxInPng(i, j)+k];
            }
        }


    //model->gpu.hssoa.RemoveDisabledAndSort(model->prms.GridY);
    stbi_image_free(png_data);
    stbi_image_free(png_data_modelled_region);

    // particle volume and mass
    model->prms.ComputeHelperVariables();

    // allocate GPU partitions
    model->gpu.initialize();
    model->gpu.split_hssoa_into_partitions();
    model->gpu.allocate_arrays();
    model->gpu.transfer_to_device();

    model->Prepare();
    spdlog::info("PreparePointsAndSetupGrid done\n");
}





void icy::SnapshotManager::ReadSnapshot(std::string fileName)
{
    /*
    {
        // Convert to std::filesystem::path for manipulation
        std::filesystem::path snPath(fileName);
        // Navigate to the target file by appending the relative path
        std::filesystem::path gridAndWindPath = snPath.parent_path() / "_grid_and_wind.h5";
        std::string grid_and_wid_Path = gridAndWindPath.string();

        H5::H5File file(grid_and_wid_Path, H5F_ACC_RDONLY);

        H5::DataSet gridDataset = file.openDataSet("grid");
        model->prms.ReadParametersFromHDF5Attributes(gridDataset);

        // Get the dataspace of the dataset
        H5::DataSpace dataspace = gridDataset.getSpace();

        // Get the dimensions of the dataset
        hsize_t dims[2];
        dataspace.getSimpleExtentDims(dims, nullptr);
        size_t gridXTotal = dims[0];
        size_t gridY = dims[1];

        if(gridXTotal != model->prms.GridXTotal || gridY != model->prms.GridY)
            throw std::runtime_error("LoadHDF5Frame grid size mismatch");

        model->gpu.hssoa.Allocate(model->prms.nPtsInitial, gridXTotal*gridY);

        // Read the dataset into the vector
        gridDataset.read(model->gpu.hssoa.grid_status_buffer.data(), H5::PredType::NATIVE_UINT8);


        // point colors

        H5::DataSet ptcDataset = file.openDataSet("point_colors_rgb");
        ptcDataset.getSpace().getSimpleExtentDims(dims, nullptr);
        if(dims[0] != model->prms.nPtsInitial) throw std::runtime_error("restore snapshot: pts count mismatch");
        ptcDataset.read(model->gpu.hssoa.point_colors_rgb.data(), H5::PredType::NATIVE_UINT32);

        model->wind_interpolator.ReadFromOwnHDF5(file);
    }

    // read state of the points
    H5::H5File file(fileName, H5F_ACC_RDONLY);
    H5::DataSet ds = file.openDataSet("pts_data");

    ds.openAttribute("SimulationStep").read(H5::PredType::NATIVE_INT, &model->prms.SimulationStep);
    ds.openAttribute("SimulationTime").read(H5::PredType::NATIVE_DOUBLE, &model->prms.SimulationTime);
    ds.openAttribute("HSSOA_size").read(H5::PredType::NATIVE_UINT, &model->gpu.hssoa.size);

    int nPtsArrays;
    ds.openAttribute("nPtsArrays").read(H5::PredType::NATIVE_INT, &nPtsArrays);

    // Get the dataspace of the dataset
    H5::DataSpace dsp = ds.getSpace();
    hsize_t dims[2];
    dsp.getSimpleExtentDims(dims, nullptr);

    if(dims[0] != SimParams::nPtsArrays || dims[1] != model->gpu.hssoa.size || nPtsArrays != dims[0])
        throw std::runtime_error("icy::SnapshotManager::ReadSnapshot array size mismatch");

    if(model->gpu.hssoa.capacity < dims[2])
        throw std::runtime_error("somehow there is not enough space in hssoa");

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

    // GPU

    model->gpu.hssoa.RemoveDisabledAndSort(model->prms.GridY); // this is actually for SpecialSelector2D

    model->gpu.initialize();
    model->gpu.split_hssoa_into_partitions();
    model->gpu.allocate_arrays();
    model->gpu.transfer_to_device();

    model->Prepare();

    spdlog::info("icy::SnapshotManager::ReadSnapshot done");
*/
}


void icy::SnapshotManager::SaveFrame(int SimulationStep, double SimulationTime)
{
    /*
    const int gridSize = model->prms.GridXTotal*model->prms.GridY;
    tmp.assign(gridSize*3,0); // in case we render a 3-channel image
    vis_vx.assign(gridSize,0);
    vis_vy.assign(gridSize,0);
    vis_r.assign(gridSize,0);
    vis_g.assign(gridSize,0);
    vis_b.assign(gridSize,0);

    vis_Jpinv.assign(gridSize,0);
    vis_P.assign(gridSize,0);
    vis_Q.assign(gridSize,0);

    vis_alpha.assign(gridSize,0);

    const int nPts = model->gpu.hssoa.size;

    for(int i=0;i<nPts;i++)
    {
        SOAIterator s = model->gpu.hssoa.begin()+i;
        //PointVector2r pos = s->getPos(model->prms.cellsize);
        int cellIdx = s->getCellIndex(model->prms.GridY);

        tmp[cellIdx]++;
        vis_vx[cellIdx] += static_cast<float>(s->getValue(SimParams::velx+0));
        vis_vy[cellIdx] += static_cast<float>(s->getValue(SimParams::velx+1));

        int pt_idx = s->getValueInt(SimParams::integer_point_idx);
        uint32_t rgb = model->gpu.hssoa.point_colors_rgb[pt_idx];
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
        int n = tmp[i];
        float coeff = 1./n;
        vis_vx[i] *= coeff;
        vis_vy[i] *= coeff;
        vis_r[i] *= coeff;
        vis_g[i] *= coeff;
        vis_b[i] *= coeff;
        vis_Jpinv[i] *= coeff;
        vis_P[i] *= coeff;
        vis_Q[i] *= coeff;
        vis_alpha[i] = std::min(n*0.2f,1.f);
    }


    std::string s_path(frame_path);
    std::filesystem::path directory_path(frame_path);
    if (!std::filesystem::exists(directory_path)) std::filesystem::create_directories(directory_path);
    int frame = SimulationStep / model->prms.UpdateEveryNthStep;

    std::string fileName = fmt::format("{}/frame_{:05d}.h5", s_path, frame);
    H5::H5File file(fileName, H5F_ACC_TRUNC);

    hsize_t dims_grid[2] = {(hsize_t)model->prms.GridXTotal, (hsize_t)model->prms.GridY};
    H5::DataSpace dataspace_grid(2, dims_grid);

    H5::DSetCreatPropList proplist;
    hsize_t chunk_dims[2] = {256, (hsize_t)model->prms.GridY};
    proplist.setChunk(2, chunk_dims);
    proplist.setDeflate(3);

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

    H5::DataSet ds_vis_r = file.createDataSet("vis_r", H5::PredType::NATIVE_FLOAT, dataspace_grid, proplist);
    ds_vis_r.write(vis_r.data(), H5::PredType::NATIVE_FLOAT);
    H5::DataSet ds_vis_g = file.createDataSet("vis_g", H5::PredType::NATIVE_FLOAT, dataspace_grid, proplist);
    ds_vis_g.write(vis_g.data(), H5::PredType::NATIVE_FLOAT);
    H5::DataSet ds_vis_b = file.createDataSet("vis_b", H5::PredType::NATIVE_FLOAT, dataspace_grid, proplist);
    ds_vis_b.write(vis_b.data(), H5::PredType::NATIVE_FLOAT);
    H5::DataSet ds_vis_alpha = file.createDataSet("vis_alpha", H5::PredType::NATIVE_FLOAT, dataspace_grid, proplist);
    ds_vis_alpha.write(vis_alpha.data(), H5::PredType::NATIVE_FLOAT);

    H5::DataSet ds_count = file.createDataSet("count", H5::PredType::NATIVE_UINT8, dataspace_grid, proplist);
    ds_count.write(tmp.data(), H5::PredType::NATIVE_UINT8);

    H5::DataSpace att_dspace(H5S_SCALAR);
    ds_count.createAttribute("SimulationStep", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &SimulationStep);
    ds_count.createAttribute("SimulationTime", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &SimulationTime);
*/
}


void icy::SnapshotManager::SaveSnapshot(int SimulationStep, double SimulationTime)
{
    /*
    //snapshot_path
    std::string s_path(snapshot_path);
    std::filesystem::path directory_path(snapshot_path);
    if (!std::filesystem::exists(directory_path)) std::filesystem::create_directories(directory_path);

    int frame = SimulationStep / model->prms.UpdateEveryNthStep;
    if(frame == 0)
    {
        // save wind data, grid/land, and point colors
        std::string fileName = s_path + "/_grid_and_wind.h5";
        H5::H5File file(fileName, H5F_ACC_TRUNC);

        // save grid/land data
        hsize_t dims_grid[2] = {(hsize_t)model->prms.GridXTotal, (hsize_t)model->prms.GridY};
        H5::DataSpace dataspace_grid(2, dims_grid);

        H5::DataSet dataset = file.createDataSet("grid", H5::PredType::NATIVE_UINT8, dataspace_grid);
        dataset.write(model->gpu.hssoa.grid_status_buffer.data(), H5::PredType::NATIVE_UINT8);


        hsize_t dims_points[1] = {model->gpu.hssoa.point_colors_rgb.size()};
        H5::DataSpace dataspace_points(1, dims_points);
        H5::DataSet dataset_pt_colors = file.createDataSet("point_colors_rgb", H5::PredType::NATIVE_UINT32, dataspace_points);
        dataset_pt_colors.write(model->gpu.hssoa.point_colors_rgb.data(), H5::PredType::NATIVE_UINT32);

        // save parameters and wind data
        model->prms.SaveParametersAsHDF5Attributes(dataset);
        model->wind_interpolator.SaveToOwnHDF5(file);
    }

    // save current state
    std::string fileName = fmt::format("{}/snapshot_{:05d}.h5", s_path, frame);
    H5::H5File file(fileName, H5F_ACC_TRUNC);


    // points
    hsize_t dims_points[2] = {SimParams::nPtsArrays, model->gpu.hssoa.size};
    H5::DataSpace dataspace_points(2, dims_points);

    // Define memory dataspace

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
*/
}






std::string icy::SnapshotManager::prepare_file_name(int gx, int gy)
{
    return fmt::format("{}/point_cache_{:05d}_{:05d}.h5", pts_cache_path, gx, gy);
}


std::pair<int, float> icy::SnapshotManager::categorizeColor(const Eigen::Vector3f& rgb)
{
//    if(rgb.x() > 0.95 && rgb.y() < 0.05 && rgb.z() < 0.05) return {-1, 0.f}; // land
//    if(rgb.x() >= 0xca/255. && rgb.y() >= 0xca/255. && rgb.z() >= 0xca/255.) return {2, 0.f}; // intact ice
//    if(rgb.x() >= 0xcd/255. && rgb.y() >= 0xd5/255. && rgb.z() >= 0xd6/255.) return {2, 0.f}; // intact ice



    // check if open water
    float minDist = std::numeric_limits<float>::max();
    int bestInterval = -1;
    float bestPosition = 0.0f;
    constexpr float openWaterThreshold = 0.05f;

    for (int i = 0; i < colordata_OpenWater.size() - 1; ++i)
    {
        // Convert std::array<float, 3> to Eigen::Vector3f
        Eigen::Vector3f p0 = arrayToEigen(colordata_OpenWater[i]);
        Eigen::Vector3f p1 = arrayToEigen(colordata_OpenWater[i + 1]);

        Eigen::Vector3f diff = p1 - p0;
        float segmentLengthSq = diff.squaredNorm();

        if (segmentLengthSq == 0.0f) continue; // Skip degenerate segments

        Eigen::Vector3f v = rgb - p0;
        Eigen::Vector3f proj = p0 + diff * v.dot(diff)/diff.squaredNorm();
        float dist = (rgb - proj).norm();

        if (dist < minDist) {
            minDist = dist;
            bestInterval = i;
        }
    }
    if(minDist < openWaterThreshold) return {0, 0.f};

    // check if crushed
    float result = projectPointOntoCurve(rgb, colordata_Solid);

    return {2,result};
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



void icy::SnapshotManager::LoadWindData(std::string fileName)
{
    /*
    spdlog::info("icy::SnapshotManager::LoadWindData - {}", fileName);
    model->prms.use_GFS_wind = true;
    model->wind_interpolator.LoadWindData(fileName,
                                          model->prms.LatMin,model->prms.LatMax,model->prms.LonMin,model->prms.LonMax,
                                          model->prms.SimulationStartUnixTime);
    model->prms.gridLatMin = model->wind_interpolator.gridLatMin;
    model->prms.gridLonMin = model->wind_interpolator.gridLonMin;
*/
}



// Function to project a point onto the curve
float icy::SnapshotManager::projectPointOntoCurve(const Eigen::Vector3f& rgb, const std::array<std::array<float, 3>, 3>& curve)
{
    float minDistance = std::numeric_limits<float>::max();
    float relativePosition = 0.0f;

    // Total length of the curve
    float totalLength = 0.0f;

    // Loop through the segments
    for (size_t i = 0; i < curve.size() - 1; ++i) {
        // Define segment points
        Eigen::Vector3f p0(curve[i][0], curve[i][1], curve[i][2]);
        Eigen::Vector3f p1(curve[i+1][0], curve[i+1][1], curve[i+1][2]);

        // Vector along the segment
        Eigen::Vector3f segment = p1 - p0;

        // Project rgb onto the infinite line
        float t = (rgb - p0).dot(segment) / segment.squaredNorm();

        // Clamp t to [0, 1] to stay on the segment
        t = std::clamp(t, 0.0f, 1.0f);

        // Compute the projected point
        Eigen::Vector3f projection = p0 + t * segment;

        // Compute the distance from rgb to the projection
        float distance = (rgb - projection).norm();

        // Update the closest segment if necessary
        if (distance < minDistance) {
            minDistance = distance;

            // Compute relative position along the curve
            float segmentLength = segment.norm();
            float positionOnSegment = totalLength + t * segmentLength;
            totalLength += segmentLength;

            relativePosition = positionOnSegment / totalLength;
        } else {
            totalLength += segment.norm();
        }
    }

    return relativePosition; // Relative position in [0, 1]
}



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

    spdlog::info("generate_and_save: attempting to generate {} points in {}x{} grid", (points_per_cell*gx*gy), gx, gy);
    buffer = thinks::PoissonDiskSampling(kRadius, kXMin, kXMax);

    const float result_ppc = (float)buffer.size()/(gx*gy);
    spdlog::info("grid: {}x{}; generated pts {:>8}; density {:>7.4}",gx, gy, buffer.size(), result_ppc);

    const float scale = sqrt(result_ppc/(points_per_cell*1.0005));
    if(scale<1.0)
    {
        spdlog::critical("requested ppc {}; generated ppc {}", points_per_cell, result_ppc);
        throw std::runtime_error("point generation error");
    }

    spdlog::info("requested ppc {}; generated ppc {}; scale {}%", points_per_cell, result_ppc, 100*(scale-1.f));

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
    spdlog::info("grid: {}x{}; updated pts {:>8}; updated_density {:>7.4}",gx, gy, buffer.size(), final_ppc);

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

    spdlog::info("sorting started");
    std::sort(buffer.begin(),buffer.end(),
              [&](std::array<float,2> &pt1, std::array<float,2> &pt2)
              {return computeIdx(pt1) < computeIdx(pt2);}
              );

    spdlog::info("sorting finished; first point: {}, {}", buffer.front()[0], buffer.front()[1]);

    std::vector<int> grid_count(gx*gy,0);
    for(auto &pt : buffer) grid_count[computeIdx(pt)]++;

    auto [it_min, it_max] = std::minmax_element(grid_count.begin(),grid_count.end());

    spdlog::info("grid fill min/max {} to {} ", *it_min, *it_max);

    std::vector<int> histogram(std::max((int)points_per_cell*3, *it_max), 0);
    for(int &count : grid_count) histogram[count]++;
    for(auto val : histogram) std::cout << val << ' ';
    std::cout << std::endl;

    // write file
    std::filesystem::path directory_path(pts_cache_path);
    if (!std::filesystem::exists(directory_path)) std::filesystem::create_directories(directory_path);
    std::string fileName = prepare_file_name(gx, gy);
    spdlog::info("saving file {}", fileName);
    H5::H5File file(fileName, H5F_ACC_TRUNC);

    const size_t nPts = buffer.size();
    hsize_t dims_pts[2] = {nPts,2};
    H5::DataSpace dataspace_pts(2, dims_pts);

    H5::DataSet dataset = file.createDataSet("coords", H5::PredType::NATIVE_FLOAT, dataspace_pts);
    dataset.write(buffer.data(), H5::PredType::NATIVE_FLOAT);

    hsize_t att_dim = 1;
    H5::DataSpace att_dspace(1, &att_dim);
    dataset.createAttribute("gx", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &gx);
    dataset.createAttribute("gy", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &gy);
    file.close();
}


bool icy::SnapshotManager::attempt_to_fill_from_cache(int gx, int gy, std::vector<std::array<float, 2>> &buffer)
{
    spdlog::info("attempting to load from cache");
    std::string fileName = prepare_file_name(gx, gy);
    std::filesystem::path file_path(fileName);
    if (!std::filesystem::exists(file_path))
    {
        spdlog::info("cached file does not exist");
        return false;
    }

    H5::H5File file(fileName, H5F_ACC_RDONLY);
    H5::DataSet dataset = file.openDataSet("coords");

    hsize_t dims[2];
    dataset.getSpace().getSimpleExtentDims(dims,NULL);
    if(dims[1] != 2) {
        spdlog::critical("something is wrong with cache file - dims[1] is {}",dims[1]);
        return false;
    }

    int saved_gx, saved_gy;
    dataset.openAttribute("gx").read(H5::PredType::NATIVE_INT, &saved_gx);
    dataset.openAttribute("gy").read(H5::PredType::NATIVE_INT, &saved_gy);
    if(saved_gx != gx || saved_gy != gy)
    {
        spdlog::critical("cache error - expected grid {}x{}; cached grid {}x{}",gx, gy, saved_gx, saved_gy);
        return false;
    }

    buffer.resize(dims[0]);
    spdlog::info("reading {} points from file", dims[0]);
    dataset.read(buffer.data(), H5::PredType::NATIVE_FLOAT);
    spdlog::info("attempt_to_fill_from_cache: finished reading from file; {} points", buffer.size());
    spdlog::info("attempt_to_fill_from_cache: first point: {}, {}", buffer.front()[0], buffer.front()[1]);

    return true;
}

