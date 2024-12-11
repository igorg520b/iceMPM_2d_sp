#include "snapshotmanager.h"
#include "model.h"

#include <spdlog/spdlog.h>
#include <H5Cpp.h>
#include <filesystem>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <utility>


void icy::SnapshotManager::LoadRawPoints(std::string fileName)
{
    spdlog::info("reading raw points file {}",fileName);
    if(!std::filesystem::exists(fileName)) throw std::runtime_error("error reading raw points file - no file");;

    H5::H5File file(fileName, H5F_ACC_RDONLY);

    H5::DataSet dataset_grid = file.openDataSet("GridLand");

    // read grid dimensions from the file
    dataset_grid.openAttribute("GridX").read(H5::PredType::NATIVE_INT, &model->prms.GridXTotal);
    dataset_grid.openAttribute("GridY").read(H5::PredType::NATIVE_INT, &model->prms.GridY);

    // compare grid_total with the size of the dataset
    hsize_t nGridNodes;
    dataset_grid.getSpace().getSimpleExtentDims(&nGridNodes, NULL);
    const int grid_total = model->prms.GridXTotal * model->prms.GridY;
    if(nGridNodes != grid_total)
    {
        spdlog::critical("nGridNodes {};  grid_total {};  grid [{} x {}]",
                         nGridNodes, grid_total, model->prms.GridXTotal, model->prms.GridY);
        throw std::runtime_error("grid size does not match the dataset in H5 file");
    }

    // read cellsize
    float cellsize;
    dataset_grid.openAttribute("cellsize").read(H5::PredType::NATIVE_FLOAT, &cellsize);
    model->prms.cellsize = (t_PointReal)cellsize;
    model->prms.cellsize_inv = 1./model->prms.cellsize;

    // get number of points
    H5::DataSet dataset_x = file.openDataSet("x");
    hsize_t nPoints;
    dataset_x.getSpace().getSimpleExtentDims(&nPoints, NULL);
    model->prms.nPtsInitial = nPoints;
    spdlog::info("nPoints {}; grid [{} x {}]", nPoints, model->prms.GridXTotal, model->prms.GridY);

    // allocate space host-side
    model->gpu.hssoa.Allocate(nPoints, grid_total);
    model->gpu.hssoa.size = nPoints;
    t_PointReal *hssoa_ptr;
    {
    std::vector<float> buffer(nPoints);

    // read
    file.openDataSet("x").read(buffer.data(), H5::PredType::NATIVE_FLOAT);
    hssoa_ptr = model->gpu.hssoa.getPointerToPosX();
    for(int i=0;i<nPoints;i++) hssoa_ptr[i] = (t_PointReal)buffer[i];

    file.openDataSet("y").read(buffer.data(), H5::PredType::NATIVE_FLOAT);
    hssoa_ptr = model->gpu.hssoa.getPointerToPosY();
    for(int i=0;i<nPoints;i++) hssoa_ptr[i] = (t_PointReal)buffer[i];

    file.openDataSet("r").read(buffer.data(), H5::PredType::NATIVE_FLOAT);
    hssoa_ptr = model->gpu.hssoa.getPointerToLine(SimParams::idx_rgb + 0);
    for(int i=0;i<nPoints;i++) hssoa_ptr[i] = (t_PointReal)buffer[i];

    file.openDataSet("g").read(buffer.data(), H5::PredType::NATIVE_FLOAT);
    hssoa_ptr = model->gpu.hssoa.getPointerToLine(SimParams::idx_rgb + 1);
    for(int i=0;i<nPoints;i++) hssoa_ptr[i] = (t_PointReal)buffer[i];

    file.openDataSet("b").read(buffer.data(), H5::PredType::NATIVE_FLOAT);
    hssoa_ptr = model->gpu.hssoa.getPointerToLine(SimParams::idx_rgb + 2);
    for(int i=0;i<nPoints;i++) hssoa_ptr[i] = (t_PointReal)buffer[i];
    }


//    {
//    std::vector<uint32_t> buffer2(nPoints);
//    dataset_grains.read(buffer2.data(), H5::PredType::NATIVE_UINT32);
//    hssoa_ptr = model->gpu.hssoa.getPointerToLine(SimParams::idx_utility_data);
//    for(int i=0;i<nPoints;i++) *reinterpret_cast<uint32_t*>(&hssoa_ptr[i]) = buffer2[i];
//    }

    // read volume attribute
    H5::Attribute att_volume = dataset_x.openAttribute("volume");
    float volume;
    att_volume.read(H5::PredType::NATIVE_FLOAT, &volume);
    model->prms.ParticleVolume = volume;

    // GRID
    dataset_grid.read(model->gpu.hssoa.grid_status_buffer.data(), H5::PredType::NATIVE_UINT8);

    file.close();


    HostSideSOA &hssoa = model->gpu.hssoa;
    hssoa.convertToIntegerCellFormat(model->prms.cellsize);

    // analyze color data
    std::vector<int> categories(colordata.size());
    for(SOAIterator it = hssoa.begin(); it!=hssoa.end(); ++it)
    {
        ProxyPoint &p = *it;
        t_PointReal r = p.getValue(SimParams::idx_rgb+0);
        t_PointReal g = p.getValue(SimParams::idx_rgb+1);
        t_PointReal b = p.getValue(SimParams::idx_rgb+2);
        Eigen::Vector3f rgb(r,g,b);
        auto [category, interpValue] = categorizeColor(rgb);
        categories[category]++;
        uint32_t val = p.getValueInt(SimParams::idx_utility_data);
        if(category <= 4)
        {
            val |= 0x20000;     // (mark for removal)
        }
        else if(category <= 11)
        {
            val |= 0x10000;     // crushed
        }
        else if(category <= 12)
        {
            val |= 0x40000;     // crushed
        }
        p.setValueInt(SimParams::idx_utility_data, val);
    }
    for(int i=0;i<categories.size();i++)
        spdlog::info("category {}; count {}", i, categories[i]);


    model->gpu.hssoa.RemoveDisabledAndSort(model->prms.GridY);
    model->gpu.hssoa.InitializeBlock();

    // particle volume and mass
    model->prms.ComputeHelperVariables();

    // allocate GPU partitions
    model->gpu.initialize();
    model->gpu.split_hssoa_into_partitions();
    model->gpu.allocate_arrays();
    model->gpu.transfer_to_device();

    model->Prepare();

    spdlog::info("LoadRawPoints done; nPoitns {}\n",nPoints);
}

std::pair<int, float> icy::SnapshotManager::categorizeColor(const Eigen::Vector3f& rgb)
{
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
