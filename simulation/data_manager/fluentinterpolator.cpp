#include "fluentinterpolator.h"

#include <iostream>
#include <filesystem>
#include <regex>
#include <vector>
#include <algorithm>
#include <spdlog/spdlog.h>
#include <fmt/core.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkVersion.h>

#include <H5Cpp.h>

namespace fs = std::filesystem;

float* FluentInterpolator::getFramePtr(int frame, int component)
{
    if(currentFrame < 0 || circularBufferIdx < 0)
        throw std::runtime_error("FluentInterpolator::getFramePtr");
    const int &gx = prms->GridXTotal;
    const int &gy = prms->GridYTotal;
    int slot = (circularBufferIdx+frame)%preloadedFrames;
    return _data.data() + gx*gy*(2*slot+component);
}

void FluentInterpolator::LoadFrame(int frame, int slot)
{
    H5::H5File file(cachedFileName, H5F_ACC_RDONLY);
    H5::DataSet ds = file.openDataSet("vx_vy");

    const int &gx = prms->GridXTotal;
    const int &gy = prms->GridYTotal;
    float *ptr = _data.data() + gx*gy*2*slot;

    hsize_t dims_mem[3] = {2, (hsize_t)gx, (hsize_t)gy};
    H5::DataSpace memspace(3, dims_mem);

    hsize_t offset[3] = {(hsize_t)frame*2, 0, 0};
    hsize_t count[3] = {2, (hsize_t)gx, (hsize_t)gy};

    H5::DataSpace dsp = ds.getSpace();
    dsp.selectHyperslab(H5S_SELECT_SET, count, offset);

    ds.read(ptr, H5::PredType::NATIVE_FLOAT, memspace, dsp);
}

void FluentInterpolator::PrepareFromCachedFile(std::string _cachedFile)
{
    this->cachedFileName = _cachedFile;
    currentFrame = -1;
    circularBufferIdx = -1;
    const int &gx = prms->GridXTotal;
    const int &gy = prms->GridYTotal;
    _data.resize(gx*gy*2*preloadedFrames); // storage to hold preloaded frames

    // check if cached data is available
    {
        H5::H5File file(cachedFileName, H5F_ACC_RDONLY);
        H5::DataSet ds = file.openDataSet("vx_vy");
        ds.openAttribute("file_count").read(H5::PredType::NATIVE_INT, &file_count);
        ds.openAttribute("interval_size").read(H5::PredType::NATIVE_INT, &interval_size);
    }
    spdlog::info("invokign SetTime(0) from PrepareFlowDataCache");
    SetTime(0);
    is_initialized = true;
    return;
}


void FluentInterpolator::PrepareFlowDataCache(std::string fileName)
{
    fs::path dirPath = fs::path(fileName).parent_path();
    std::string baseName = fs::path(fileName).stem().string(); // Extract "fileName.cas"
    baseName = baseName.substr(0, baseName.size() - 4); // Remove ".cas"
    geometryFilePrefix = dirPath.string() + "/" + baseName;

    currentFrame = -1;
    circularBufferIdx = -1;

    const int &gx = prms->GridXTotal;
    const int &gy = prms->GridYTotal;
    _data.resize(gx*gy*2*preloadedFrames); // storage to hold preloaded frames

    // check if cached data is available
    cachedFileName = CachedFileName();
    if(fs::exists(cachedFileName))
    {
        {
            H5::H5File file(cachedFileName, H5F_ACC_RDONLY);
            H5::DataSet ds = file.openDataSet("vx_vy");
            ds.openAttribute("file_count").read(H5::PredType::NATIVE_INT, &file_count);
            ds.openAttribute("interval_size").read(H5::PredType::NATIVE_INT, &interval_size);
        }
        spdlog::info("invokign SetTime(0) from PrepareFlowDataCache");
        SetTime(0);
        is_initialized = true;
        return;
    }
    spdlog::info("preparing cache file {}", cachedFileName);

    // generate and save HDF5 file with rasterized flow data

    std::regex filePattern(baseName + "-(\\d{5})\\.dat\\.h5");

    std::vector<int> numbers;   // file/frame numbers

    for (const auto& entry : fs::directory_iterator(dirPath)) {
        std::smatch match;
        std::string fName = entry.path().filename().string();
        if (std::regex_match(fName, match, filePattern)) {
            numbers.push_back(std::stoi(match[1].str()));
        }
    }

    if (numbers.empty()) {
        interval_size = 0;
        file_count = 0;
        return;
    }

    std::sort(numbers.begin(), numbers.end());
    interval_size = (numbers.size() > 1) ? (numbers[1] - numbers[0]) : 0;
    file_count = static_cast<int>(numbers.size())-1;


    // prepare the cache file
    H5::H5File file(cachedFileName, H5F_ACC_TRUNC);

    hsize_t dims[3] = {(hsize_t)file_count*2, (hsize_t)gx, (hsize_t)gy};
    H5::DataSpace dataspace(3, dims);
    H5::DataSet dataset = file.createDataSet("vx_vy", H5::PredType::NATIVE_FLOAT, dataspace);

    H5::DataSpace att_dspace(H5S_SCALAR);
    dataset.createAttribute("gx", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &gx);
    dataset.createAttribute("gy", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &gy);
    dataset.createAttribute("file_count", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &file_count);
    dataset.createAttribute("interval_size", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &interval_size);

    hsize_t dims_mem[3] = {2, (hsize_t)gx, (hsize_t)gy};
    H5::DataSpace memspace(3, dims_mem);

    // rasterize frames one by one into HDF5
    for(int i=0;i<file_count;i++)
    {
        size_t frameSize = gx*gy;
        spdlog::info("FluentInterpolator::ScanDirectory: loading frame {}",i);
//        ParseDataFrame(i+1, _data.data(), _data.data()+frameSize);
        ParseDataFrame(file_count/2, _data.data(), _data.data()+frameSize);

        // define hyperslab
        hsize_t offset[3] = {(hsize_t)i*2, 0, 0};
        hsize_t count[3] = {2, (hsize_t)gx, (hsize_t)gy};
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);
        dataset.write(_data.data(),H5::PredType::NATIVE_FLOAT,memspace,dataspace);
    }

    spdlog::info("FluentInterpolator::ScanDirectory: files {}; interval {}; prefix {}",
                 file_count, interval_size, geometryFilePrefix);
    is_initialized = true;

    SetTime(0);
}




void FluentInterpolator::ParseDataFrame(int frame, float *U, float *V)
{
    const int &gx = prms->GridXTotal;
    const int &gy = prms->GridYTotal;

    const double &scale = prms->FluentDataScale;
    const double &offsetX = prms->FluentDataOffsetX;
    const double &offsetY = prms->FluentDataOffsetY;

    const int &ox = prms->ModeledRegionOffsetX;
    const int &oy = prms->ModeledRegionOffsetY;

    const int &vgx = prms->InitializationImageSizeX;
    const int &vgy = prms->InitializationImageSizeY;

    const double &h = prms->cellsize;

    // FFF.1-3-00130.dat.h5
    // FFF.1-3.cas.h5
    //    std::string dataFile1 = fmt::format("{}-{:05}.dat.h5", geometryFilePrefix, frame*interval);
    std::string dataFileName = fmt::format("{}-{:05}.dat.h5", geometryFilePrefix, frame*interval_size);
    std::string casFileName = fmt::format("{}.cas.h5", geometryFilePrefix);
    spdlog::info("data file {}", dataFileName);
    spdlog::info("cas file {}", casFileName);



    vtkNew<vtkFLUENTCFFCustomReader> fluentReader;
    vtkNew<vtkTransform> transform;
    vtkNew<vtkTransformFilter> transformFilter;
    vtkNew<vtkCellDataToPointData> filter_cd2pd;
    vtkNew<vtkImageData> imageData;
    vtkNew<vtkProbeFilter> probeFilter;



    fluentReader->SetDataFileName(dataFileName.c_str());  // Set the mesh file (.cas.h5)
    fluentReader->SetFileName(casFileName.c_str());  // Set the mesh file (.cas.h5)
    fluentReader->Update(); // Update the reader

    vtkUnstructuredGrid* extractedGrid = vtkUnstructuredGrid::SafeDownCast(fluentReader->GetOutput()->GetBlock(0));
//    extractedGrid->GetCellData()->SetActiveScalars("SV_U");

    transform->Identity();
    transform->Scale(scale, scale, scale);  // Scale uniformly
    transform->Translate(offsetX, offsetY, 0.);  // Offset

    // Apply the transform using vtkTransformFilter
    transformFilter->SetTransform(transform);
    transformFilter->SetInputData(extractedGrid);
    transformFilter->Update();

    filter_cd2pd->SetInputData(transformFilter->GetUnstructuredGridOutput());
    filter_cd2pd->Update();

    // set the size of vtkImageData
    imageData->SetDimensions(vgx, vgy, 1);
    imageData->SetSpacing(h, h, 1.0); // Adjust spacing as needed
//    imageData->SetOrigin(bounds[0], bounds[2], 0.0); // Adjust origin as needed

    // transfer to imageData via probeFilter
    probeFilter->SetInputData(imageData); // Set the 2D structured grid as input
    probeFilter->SetSourceData(filter_cd2pd->GetOutput()); // Set the unstructured grid as source
    probeFilter->Update();

    vtkImageData* probedData = vtkImageData::SafeDownCast(probeFilter->GetOutput());

    vtkDoubleArray* dataArray = vtkDoubleArray::SafeDownCast(probedData->GetPointData()->GetScalars("SV_V"));


    if (!dataArray) {
        spdlog::error("Failed to extract SV_V from probedData");
        throw std::runtime_error("Failed to extract SV_V from probedData");
    }
    for(int i=0;i<gx;i++)
        for(int j=0;j<gy;j++)
        {
            int idx1 = j + gy*i;
            int idx2 = (ox+i) + vgx * (oy+j);
            V[idx1] = dataArray->GetValue(idx2);
            if(gpu->grid_status_buffer[idx1]==0) V[idx1]=0;
        }

    dataArray = vtkDoubleArray::SafeDownCast(probedData->GetPointData()->GetScalars("SV_U"));
    if (!dataArray) {
        spdlog::error("Failed to extract SV_U from probedData");
        throw std::runtime_error("Failed to extract SV_U from probedData");
    }
    for(int i=0;i<gx;i++)
        for(int j=0;j<gy;j++)
        {
            int idx1 = j + gy*i;
            int idx2 = (ox+i) + vgx * (oy+j);
            U[idx1] = dataArray->GetValue(idx2);
            if(gpu->grid_status_buffer[idx1]==0) U[idx1]=0;
        }

    // smooth the data
    applyDiffusion(V, gx, gy, 0.1, 1.0, 50);
    applyDiffusion(U, gx, gy, 0.1, 1.0, 50);
}


Eigen::Vector2f FluentInterpolator::getInterpolation(int i, int j) const
{
    const int &gx = prms->GridXTotal;
    const int &gy = prms->GridYTotal;

    const int idx = j + gy*i;
    const float *fd1u = _data.data()+(circularBufferIdx*2+0)*gx*gy;
    const float *fd1v = _data.data()+(circularBufferIdx*2+1)*gx*gy;
    const float *fd2u = _data.data()+(((circularBufferIdx+1)%preloadedFrames)*2+0)*gx*gy;
    const float *fd2v = _data.data()+(((circularBufferIdx+1)%preloadedFrames)*2+1)*gx*gy;

    const float vx = (1-position)*fd1u[idx] + position*fd2u[idx];
    const float vy = (1-position)*fd1v[idx] + position*fd2v[idx];
    return Eigen::Vector2f(vx, vy);
}


bool FluentInterpolator::SetTime(double t)
{
    //    spdlog::info("FluentInterpolator::SetTime {}", t);
    this->currentTime = t;
    const double frame_interval = interval_size * prms->FrameTimeInterval;

    // which file indices to read
    int rawInterval = static_cast<int>(std::floor(currentTime / frame_interval));

    // Wrap around using modulo operation

    const int frameFrom = (rawInterval % file_count);
//    const int frameTo = ((rawInterval+1) % file_count);

    bool frameChanged = false;
    if(currentFrame > 0 && frameFrom == (currentFrame+1)%file_count)
    {
        frameChanged = true;
        circularBufferIdx = (circularBufferIdx+1)%preloadedFrames;
        currentFrame = frameFrom;

        // Launch LoadFrame asynchronously, but ensure only one instance runs at a time
        {
            std::unique_lock<std::mutex> lock(loadMutex);
            loadCV.wait(lock, [this] { return !loadInProgress.load(); });
            loadInProgress = true;
        }

        if (loadThread && loadThread->joinable()) {
            loadThread->join();
        }

        loadThread = std::thread([this, frameFrom]() {
            spdlog::info("async loading of frame {}",(frameFrom + 2) % file_count);
            LoadFrame((frameFrom + 2) % file_count, (circularBufferIdx + 2) % preloadedFrames);

            {
                std::lock_guard<std::mutex> lock(loadMutex);
                loadInProgress = false;
            }
            loadCV.notify_one();
        });

//        LoadFrame((currentFrame+2)%file_count, (circularBufferIdx+2)%preloadedFrames); // needs to be invoked asynchronously
    }
    if(frameFrom != currentFrame)
    {
        frameChanged = true;
        currentFrame = frameFrom;
        circularBufferIdx = 0;
        // reload circular buffer entirely
        spdlog::info("complete reloading of frame buffer {}",currentFrame);
        for(int i=0;i<preloadedFrames;i++) LoadFrame((currentFrame+i)%file_count, i);
    }
    else
    {
        // no change in the interval
    }

    // Compute position within the interval (0 to 1)
    double intervalStartTime = rawInterval * frame_interval;
    position = (currentTime - intervalStartTime) / frame_interval; // 0 to 1 range

    return frameChanged;
}




std::string FluentInterpolator::CachedFileName()
{
    const int &gx = prms->GridXTotal;
    const int &gy = prms->GridYTotal;
    const int &ox = prms->ModeledRegionOffsetX;
    const int &oy = prms->ModeledRegionOffsetY;
    const int &vgx = prms->InitializationImageSizeX;
    const int &vgy = prms->InitializationImageSizeY;
    std::string result = fmt::format("{}/flc_{:05}_{:05}_{:05}_{:05}_{:05}_{:05}.h5", flow_cache_path, gx, gy, ox, oy, vgx, vgy);
    return result;
}

FluentInterpolator::FluentInterpolator()
{
    spdlog::info("FluentInterpolator::FluentInterpolator()");
    std::cout << "VTK Version: " << VTK_MAJOR_VERSION << "."
              << VTK_MINOR_VERSION << "."
              << VTK_BUILD_VERSION << std::endl;

    std::filesystem::path directory_path(flow_cache_path);
    if (!std::filesystem::exists(directory_path)) std::filesystem::create_directories(directory_path);
}

FluentInterpolator::~FluentInterpolator()
{
    if (loadThread && loadThread->joinable()) {
        loadThread->join();
    }
}

void FluentInterpolator::TestLoad(double scale, double ox, double oy)
{
    std::string dataFileName = fmt::format("{}-{:05}.dat.h5", geometryFilePrefix, 0);
    std::string casFileName = fmt::format("{}.cas.h5", geometryFilePrefix);

    vtkNew<vtkFLUENTCFFCustomReader> fluentReader;
    fluentReader->SetDataFileName(dataFileName.c_str());  // Set the mesh file (.cas.h5)
    fluentReader->SetFileName(casFileName.c_str());  // Set the mesh file (.cas.h5)
    fluentReader->Update(); // Update the reader

    vtkUnstructuredGrid* extractedGrid = vtkUnstructuredGrid::SafeDownCast(fluentReader->GetOutput()->GetBlock(0));

    vtkNew<vtkTransform> transform;
    vtkNew<vtkTransformFilter> transformFilter;

    transform->Identity();
    transform->Scale(scale, scale, scale);  // Scale uniformly
    transform->Translate(ox, oy, 0.);  // Offset

    // Apply the transform using vtkTransformFilter
    transformFilter->SetTransform(transform);
    transformFilter->SetInputData(extractedGrid);
    transformFilter->Update();


    mapper->SetInputData(transformFilter->GetUnstructuredGridOutput());

    actor_original->SetMapper(mapper);
    actor_original->GetProperty()->SetRepresentationToWireframe();
    actor_original->GetProperty()->LightingOff();
    actor_original->PickableOff();
    actor_original->GetProperty()->SetColor(0.7,0.1,0.1);
}

void FluentInterpolator::applyDiffusion(float* V, int gx, int gy, float D, float dt, int steps) {
    // Create a temporary array to store the updated values
    std::vector<float> temp(gx * gy, 0.0f);

    // Loop over the number of diffusion steps
    for (int step = 0; step < steps; step++) {
        // Apply diffusion to each grid point
        for (int i = 0; i < gx; i++) {
            for (int j = 0; j < gy; j++) {
                int idx = j + gy * i;

                // Get neighboring indices (with boundary checks)
                int i_prev = (i > 0) ? i - 1 : i;
                int i_next = (i < gx - 1) ? i + 1 : i;
                int j_prev = (j > 0) ? j - 1 : j;
                int j_next = (j < gy - 1) ? j + 1 : j;

                // Compute second derivatives using finite differences
                float d2u_dx2 = V[j + gy * i_prev] - 2 * V[idx] + V[j + gy * i_next];
                float d2u_dy2 = V[j_prev + gy * i] - 2 * V[idx] + V[j_next + gy * i];

                // Update the value using the diffusion equation
                temp[idx] = V[idx] + D * dt * (d2u_dx2 + d2u_dy2);
            }
        }

        // Copy updated values back to the original array
        std::copy(temp.begin(), temp.end(), V);
        for (int i = 0; i < gx; i++) {
            for (int j = 0; j < gy; j++) {
                int idx = j + gy * i;
                if(gpu->grid_status_buffer[idx]==0) V[idx]=0;
            }
        }


    }
}
