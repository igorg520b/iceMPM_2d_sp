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

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


void icy::SnapshotManager::load_png(std::string pngFileName)
{
    if (!std::filesystem::exists(pngFileName))
    {
        spdlog::critical("png file does not exist: {}", pngFileName);
        throw std::runtime_error("png does not exist");
    }

    const char* filename = pngFileName.c_str();
    png_data = stbi_load(filename, &imgx, &imgy, &channels, 3); // request 3 channels - RGB

    if (!png_data)
    {
        spdlog::critical("Failed to load image: {}", filename);
        throw std::runtime_error("png not loaded");
    }

    model->prms.GridXTotal = imgx;
    model->prms.GridY = imgy;

    spdlog::info("image ({} x {}) loaded; {} channels", imgx, imgy, channels);
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


std::string icy::SnapshotManager::prepare_file_name(int gx, int gy)
{
    std::string fileName = fmt::format("{}/point_cache_{:05d}_{:05d}.h5", pts_cache_path, gx, gy);
    return fileName;

}


std::pair<int, float> icy::SnapshotManager::categorizeColor(const Eigen::Vector3f& rgb)
{
    if(rgb.x() > 0.95 && rgb.y() < 0.05 && rgb.z() < 0.05) return {-1, 0}; // land


    float minDist = std::numeric_limits<float>::max();
    int bestInterval = -1;
    float bestPosition = 0.0f;

    for (size_t i = 0; i < colordata.size() - 1; ++i) {
        // Convert std::array<float, 3> to Eigen::Vector3f
        Eigen::Vector3f p0 = arrayToEigen(colordata[i]);
        Eigen::Vector3f p1 = arrayToEigen(colordata[i + 1]);

        Eigen::Vector3f diff = p1 - p0;
        float segmentLengthSq = diff.squaredNorm();

        if (segmentLengthSq == 0.0f) continue; // Skip degenerate segments

        Eigen::Vector3f v = rgb - p0;
        float t = v.dot(diff) / segmentLengthSq;
        t = std::clamp(t, 0.0f, 1.0f);

        Eigen::Vector3f closestPoint = p0 + t * diff;
        float distSq = (rgb - closestPoint).squaredNorm();

        if (distSq < minDist) {
            minDist = distSq;
            bestInterval = static_cast<int>(i);
            bestPosition = t;
        }
    }

    return {bestInterval, bestPosition};
}

void icy::SnapshotManager::PreparePointsAndSetupGrid(std::string fileName)
{
    spdlog::info("icy::SnapshotManager::PreparePointsAndSetupGrid {}",fileName);
    load_png(fileName);

    // either load or generate points
    std::vector<std::array<float, 2>> buffer;   // initial buffer for points
    bool cache_result = attempt_to_fill_from_cache(imgx, imgy, buffer);
    if(!cache_result) generate_and_save(imgx, imgy, 5.0, buffer);

    // function for obtaining index in png_data from the pixel's 2D index (i,j)
    auto idxInPng = [&](int px, int py) -> int { return 3*((imgy - py - 1) * imgx + px); };

    // convert unscaled point coordinates to (i,j) pair on the grid
    auto idxPt = [&](std::array<float,2> &pt) -> std::pair<int,int> {
        return {(int)((double)pt[0]*(imgx-1) + 0.5), (int)((double)pt[1]*(imgx-1) + 0.5)};
    };

    // now that we have png_data, classify points and count those remaining (not land or open water)
    model->prms.nPtsInitial = std::count_if(buffer.begin(), buffer.end(),
                                               [&](std::array<float,2> &pt) {
        auto [i,j] = idxPt(pt);
        if(i<=1 || j<=1 || i>= (imgx-2) || j>= (imgy-2)) return false; // exclude points right at the boundary
        int idx_png = idxInPng(i,j);
        unsigned char r = png_data[idx_png + 0];
        unsigned char g = png_data[idx_png + 1];
        unsigned char b = png_data[idx_png + 2];
        if(r > 250 && g < 5 && b < 5) return false;
        Eigen::Vector3f rgb((float)r/255.,(float)g/255.,(float)b/255.);
        auto [category, interpValue] = categorizeColor(rgb);
        if(category <= 4) return false; // either land or water
        else return true; // some sort of ice
    });

    // allocate hssoa + grid
    model->gpu.hssoa.Allocate(model->prms.nPtsInitial, imgx*imgy);
    model->gpu.hssoa.size = model->prms.nPtsInitial;
    spdlog::info("PreparePointsAndSetupGrid: nPoints {}; grid [{} x {}]",
                 model->prms.nPtsInitial, model->prms.GridXTotal, model->prms.GridY);

    // compute the X-dimension from lat/lon using haversineDistance()
    double XScale = haversineDistance((model->prms.LatMin + model->prms.LatMax)*0.5, model->prms.LonMin, model->prms.LonMax);

    // compute cellsize
    model->prms.cellsize = XScale / (model->prms.GridXTotal-1);
    model->prms.cellsize_inv = 1.0/model->prms.cellsize;
    spdlog::info("assumed horizontal scale: {}; cellsize {}", XScale, model->prms.cellsize);

    // sequentially transfer points from "buffer" to hssoa (only those counted), scale coordinates
    int count = 0;
    HostSideSOA &hssoa = model->gpu.hssoa;
    for(std::array<float,2> &pt : buffer)
    {

        auto [i,j] = idxPt(pt);
        if(i<=1 || j<=1 || i>= (imgx-2) || j>= (imgy-2)) continue; // exclude points at the boundary
        int idx_png = idxInPng(i,j);
        unsigned char r = png_data[idx_png + 0];
        unsigned char g = png_data[idx_png + 1];
        unsigned char b = png_data[idx_png + 2];
        Eigen::Vector3f rgb((float)r/255.,(float)g/255.,(float)b/255.);
        auto [category, interpValue] = categorizeColor(rgb);
        if(category <= 4) continue; // either land or water

        // write into SOA
        if(count == model->prms.nPtsInitial)
        {
            spdlog::critical("when transferring from buffer to SOA, there is and index error");
            throw std::runtime_error("when transferring from buffer to SOA, there is and index error");
        }
        SOAIterator it = hssoa.begin()+count;
        ProxyPoint &p = *it;
        p.setValue(SimParams::posx, pt[0]*XScale);      // set point's x-position in SOA
        p.setValue(SimParams::posx+1, pt[1]*XScale);    // set point's y-position in SOA
        p.setValue(SimParams::idx_rgb+0, rgb[0]);
        p.setValue(SimParams::idx_rgb+1, rgb[1]);
        p.setValue(SimParams::idx_rgb+2, rgb[2]);

        uint32_t val = 0;
        if(category <= 11) val |= 0x10000;     // crushed
        else if(category <= 12) val |= 0x40000;     // weakened
        p.setValueInt(SimParams::idx_utility_data, val);

        count++;
    }

    // convert points to cell-based local coordinates
    hssoa.convertToIntegerCellFormat(model->prms.cellsize);

    model->prms.ParticleVolume = XScale*XScale * (float)imgy/imgx / buffer.size();


    // transfer / set grid land bit
    for(int i=0;i<imgx;i++)
        for(int j=0;j<imgy;j++)
        {
            int idx_png = idxInPng(i,j);
            unsigned char r = png_data[idx_png + 0];
            unsigned char g = png_data[idx_png + 1];
            unsigned char b = png_data[idx_png + 2];

            bool is_land = (r > 250 && g < 5 & b < 5);
            hssoa.grid_status_buffer[j + i*imgy] = is_land ? 1 : 0;
        }

    //model->gpu.hssoa.RemoveDisabledAndSort(model->prms.GridY);
    model->gpu.hssoa.InitializeBlock();

    // particle volume and mass
    model->prms.ComputeHelperVariables();

    // allocate GPU partitions
    model->gpu.initialize();
    model->gpu.split_hssoa_into_partitions();
    model->gpu.allocate_arrays();
    model->gpu.transfer_to_device();

    model->Prepare();

    stbi_image_free(png_data);
    spdlog::info("PreparePointsAndSetupGrid done\n");
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
    spdlog::info("icy::SnapshotManager::LoadWindData - {}", fileName);

    H5::H5File file(fileName, H5F_ACC_RDONLY);
    H5::DataSet dataset_valid_time = file.openDataSet("valid_time");
    H5::DataSpace dataspace = dataset_valid_time.getSpace();

    // Get the size of the dataset
    hsize_t datasetSize;
    dataspace.getSimpleExtentDims(&datasetSize, nullptr);

    std::vector<int64_t> &valid_time = model->valid_time;
    valid_time.resize(datasetSize);
    dataset_valid_time.read(valid_time.data(), H5::PredType::NATIVE_LLONG);
    spdlog::info("valid_time read {} values", datasetSize);



    spdlog::info("icy::SnapshotManager::LoadWindData done");


    model->use_GFS_wind = true;

}
