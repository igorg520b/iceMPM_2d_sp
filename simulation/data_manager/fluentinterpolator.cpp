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


namespace fs = std::filesystem;


void FluentInterpolator::TestLoad(double scale, double ox, double oy)
{
    std::string dataFileName = R"(/home/s2/Documents/vm_shared/fluent4/data2/FFF.1-3-01100.dat.h5)";
    fluentReader->SetDataFileName(dataFileName.c_str());  // Set the mesh file (.cas.h5)

    std::string fileName = R"(/home/s2/Documents/vm_shared/fluent4/data2/FFF.1-3.cas.h5)";
    fluentReader->SetFileName(fileName.c_str());  // Set the mesh file (.cas.h5)
    fluentReader->Update(); // Update the reader

    vtkUnstructuredGrid* extractedGrid = vtkUnstructuredGrid::SafeDownCast(fluentReader->GetOutput()->GetBlock(0));
    double bounds[6];
    extractedGrid->GetBounds(bounds);
    spdlog::info("bounds x [{}, {}]; y[{}, {}]", bounds[0], bounds[1], bounds[2], bounds[3]);
    double width = bounds[1]-bounds[0];
    double height = bounds[3]-bounds[2];

    extractedGrid->GetCellData()->SetActiveScalars("SV_U");


    transform->Identity();
    transform->Scale(scale, scale, scale);  // Scale uniformly
    transform->Translate(ox, oy, 0.);  // Offset

    // Apply the transform using vtkTransformFilter
    transformFilter->SetTransform(transform);
    transformFilter->SetInputData(extractedGrid);
    transformFilter->Update();

    filter_cd2pd->SetInputData(transformFilter->GetUnstructuredGridOutput());
    filter_cd2pd->Update();

    mapper->SetInputConnection(filter_cd2pd->GetOutputPort());
//    mapper->ScalarVisibilityOn();
    mapper->ScalarVisibilityOff();
    mapper->SetLookupTable(lut);
    mapper->SetScalarRange(-0.1, 0.1);

    actor_original->SetMapper(mapper);
    actor_original->GetProperty()->SetRepresentationToWireframe();
    actor_original->GetProperty()->LightingOff();
    actor_original->PickableOff();
    actor_original->GetProperty()->SetColor(0.7,0.1,0.1);


/*
    // transfer to imageData via probeFilter
    const int gx = 1000;
    const int gy = (int)(gx*height/width);
    const double spacing = width/(gx-1);

    imageData->SetDimensions(gx, gy, 1);
    imageData->SetSpacing(spacing, spacing, 1.0); // Adjust spacing as needed
    imageData->SetOrigin(bounds[0], bounds[2], 0.0); // Adjust origin as needed

    probeFilter->SetInputData(imageData); // Set the 2D structured grid as input
    probeFilter->SetSourceData(filter_cd2pd->GetOutput()); // Set the unstructured grid as source
    probeFilter->Update();

    vtkImageData* probedData = vtkImageData::SafeDownCast(probeFilter->GetOutput());

    // transfer the image into array

    vtkDoubleArray* dataArray = vtkDoubleArray::SafeDownCast(probedData->GetPointData()->GetScalars("SV_U"));
    if (!dataArray) {
        spdlog::error("Failed to extract SV_U from probedData");
        throw std::runtime_error("Failed to extract SV_U from probedData");
    }

    flatData_U.resize(gx*gy);
    for (int idx = 0; idx < gx * gy; ++idx) flatData_U[idx] = dataArray->GetValue(idx);

    flatData_V.resize(gx*gy);
    dataArray = vtkDoubleArray::SafeDownCast(probedData->GetPointData()->GetScalars("SV_V"));
    if (!dataArray) {
        spdlog::error("Failed to extract SV_V from probedData");
        throw std::runtime_error("Failed to extract SV_V from probedData");
    }
    for (int idx = 0; idx < gx * gy; ++idx) flatData_V[idx] = dataArray->GetValue(idx);


    imageData->AllocateScalars(VTK_FLOAT,1);
    float* data = static_cast<float*>(imageData->GetScalarPointer());
    for (int idx = 0; idx < gx * gy; ++idx) data[idx] = flatData_V[idx];
    imageData->Modified();

//    probeMapper->SetInputData(probedData);
    probeMapper->SetInputData(imageData);
    probeMapper->SetScalarRange(-0.1, 0.1);
    probeMapper->ScalarVisibilityOn();
    probeMapper->SetLookupTable(lut);

    actor->SetMapper(probeMapper);


//    actor->SetMapper(mapper);
    actor->GetProperty()->SetEdgeVisibility(true);
    actor->GetProperty()->SetEdgeColor(0.2,0.2,0.2);
    actor->GetProperty()->LightingOff();
    actor->PickableOff();
    actor->GetProperty()->SetColor(0.1,0.1,0.1);



    //    vtkMultiBlockDataSet* set = fluentReader->GetOutput();
    //    spdlog::info("FluentInterpolator::LoadDataFrame; dataFile {}, blocks {}; cells {}; points {}", fileName,
    //                 set->GetNumberOfBlocks(), set->GetNumberOfCells(), set->GetNumberOfPoints());

    //    vtkCellData* cellData = extractedGrid->GetCellData();
    //    if (!cellData || cellData->GetNumberOfArrays() == 0) {
    //        std::cerr << "Error: No per-cell scalar data found!" << std::endl;
    //        throw std::runtime_error("celldata");
    //    }
    //   spdlog::info("celldata arrays {}", cellData->GetNumberOfArrays());
*/
}



FluentInterpolator::FluentInterpolator()
{
    spdlog::info("FluentInterpolator::FluentInterpolator()");
    std::cout << "VTK Version: " << VTK_MAJOR_VERSION << "."
              << VTK_MINOR_VERSION << "."
              << VTK_BUILD_VERSION << std::endl;

}


void FluentInterpolator::ScanDirectory(std::string fileName)
{
    fs::path dirPath = fs::path(fileName).parent_path();
    std::string baseName = fs::path(fileName).stem().string(); // Extract "fileName.cas"
    baseName = baseName.substr(0, baseName.size() - 4); // Remove ".cas"
    geometryFilePrefix = dirPath.string() + "/" + baseName;

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
    file_count = static_cast<int>(numbers.size());


    // load all these into memory
    const int &gx = prms->GridXTotal;
    const int &gy = prms->GridYTotal;

    _data.resize(file_count*gx*gy*2);

    for(int i=0;i<file_count;i++)
    {
        spdlog::info("FluentInterpolator::ScanDirectory: loading frame {}",i);
        LoadDataFrame(i, getData(i,0), getData(i,1));
    }



    spdlog::info("FluentInterpolator::ScanDirectory: files {}; interval {}; prefix {}",
                 file_count, interval_size, geometryFilePrefix);
    is_initialized = true;

    SetTime(0);
}




void FluentInterpolator::LoadDataFrame(int frame, float *U, float *V)
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
        }
}


Eigen::Vector2f FluentInterpolator::getInterpolation(int i, int j)
{
    const int &gx = prms->GridXTotal;
    const int &gy = prms->GridYTotal;

    int idx = j + gy*i;
    float *fd1u = getData(currentFrame, 0);
    float *fd1v = getData(currentFrame, 1);
    float *fd2u = getData(currentFrame, 2);
    float *fd2v = getData(currentFrame, 3);

    float vx = (1-position)*fd1u[idx] + position*fd2u[idx];
    float vy = (1-position)*fd1v[idx] + position*fd2v[idx];
    return Eigen::Vector2f(vx, vy);
}

float* FluentInterpolator::getData(int frame, int subIdx)
{
    const int &gx = prms->GridXTotal;
    const int &gy = prms->GridYTotal;
    const int framesize = gx*gy;
    return _data.data() + framesize*((frame*2+subIdx)%file_count);
}

bool FluentInterpolator::SetTime(double t)
{
    //    spdlog::info("FluentInterpolator::SetTime {}", t);
    this->currentTime = t;

    // which file indices to read
    int rawInterval = static_cast<int>(std::floor(currentTime / interval_size));

    // Wrap around using modulo operation

    const int frameFrom = (rawInterval % file_count);
    const int frameTo = ((rawInterval+1) % file_count);

    bool frameChanged = false;
    if(frameFrom != currentFrame)
    {
        frameChanged = true;
        currentFrame = frameFrom;
    }

    // Compute position within the interval (0 to 1)
    double intervalStartTime = rawInterval * interval_size;
    position = (currentTime - intervalStartTime) / interval_size; // 0 to 1 range

    return frameChanged;
}

