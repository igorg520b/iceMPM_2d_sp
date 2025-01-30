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


void FluentInterpolator::TestLoad()
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

    filter_cd2pd->SetInputData(extractedGrid);
    filter_cd2pd->Update();

    mapper->SetInputConnection(filter_cd2pd->GetOutputPort());
//    mapper->ScalarVisibilityOn();
    mapper->ScalarVisibilityOff();
    mapper->SetLookupTable(lut);
    mapper->SetScalarRange(-0.1, 0.1);


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


    actor_original->SetMapper(mapper);
    actor_original->GetProperty()->SetRepresentationToWireframe();
    actor_original->GetProperty()->LightingOff();
    actor_original->PickableOff();
    actor_original->GetProperty()->SetColor(0.1,0.1,0.1);

    //    vtkMultiBlockDataSet* set = fluentReader->GetOutput();
    //    spdlog::info("FluentInterpolator::LoadDataFrame; dataFile {}, blocks {}; cells {}; points {}", fileName,
    //                 set->GetNumberOfBlocks(), set->GetNumberOfCells(), set->GetNumberOfPoints());

    //    vtkCellData* cellData = extractedGrid->GetCellData();
    //    if (!cellData || cellData->GetNumberOfArrays() == 0) {
    //        std::cerr << "Error: No per-cell scalar data found!" << std::endl;
    //        throw std::runtime_error("celldata");
    //    }
    //   spdlog::info("celldata arrays {}", cellData->GetNumberOfArrays());

}



FluentInterpolator::FluentInterpolator()
{
    spdlog::info("FluentInterpolator::FluentInterpolator()");
    std::cout << "VTK Version: " << VTK_MAJOR_VERSION << "."
              << VTK_MINOR_VERSION << "."
              << VTK_BUILD_VERSION << std::endl;

}

void FluentInterpolator::LoadDataFrame(int frame)
{
    // FFF.1-3-00130.dat.h5
    // FFF.1-3.cas.h5
//    std::string dataFile1 = fmt::format("{}-{:05}.dat.h5", geometryFilePrefix, frame*interval);
    std::string dataFile1 = fmt::format("{}-{:05}", geometryFilePrefix, frame*interval);
    spdlog::info("FluentInterpolator::LoadDataFrame; opening {}", dataFile1);
/*
    fluentReader->SetFileName(dataFile1.c_str());  // Set the mesh file (.cas.h5)
    fluentReader->Update(); // Update the reader

    vtkMultiBlockDataSet* set = fluentReader->GetOutput();
    spdlog::info("FluentInterpolator::LoadDataFrame; frame {}; dataFile {}, blocks {}; cells {}; points {}", frame, dataFile1,
        set->GetNumberOfBlocks(), set->GetNumberOfCells(), set->GetNumberOfPoints());
*/
}



void FluentInterpolator::ScanDirectory(std::string fileName)
{
    fs::path dirPath = fs::path(fileName).parent_path();
    std::string baseName = fs::path(fileName).stem().string(); // Extract "fileName.cas"
    baseName = baseName.substr(0, baseName.size() - 4); // Remove ".cas"
    geometryFilePrefix = dirPath.string() + "/" + baseName;

    std::regex filePattern(baseName + "-(\\d{5})\\.dat\\.h5");
    std::vector<int> numbers;

    for (const auto& entry : fs::directory_iterator(dirPath)) {
        std::smatch match;
        std::string fName = entry.path().filename().string();
        if (std::regex_match(fName, match, filePattern)) {
            numbers.push_back(std::stoi(match[1].str()));
        }
    }

    if (numbers.empty()) {
        interval = 0;
        file_count = 0;
        return;
    }

    std::sort(numbers.begin(), numbers.end());
    interval = (numbers.size() > 1) ? (numbers[1] - numbers[0]) : 0;
    file_count = static_cast<int>(numbers.size());

    spdlog::info("FluentInterpolator::ScanDirectory: files {}; interval {}; prefix {}", file_count, interval, geometryFilePrefix);
}

