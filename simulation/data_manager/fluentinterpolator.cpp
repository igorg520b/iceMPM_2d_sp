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
    std::string fileName = R"(/home/s2/Documents/vm_shared/fluent4/data2/FFF.1-3.cas.h5)";
    fluentReader->SetFileName(fileName.c_str());  // Set the mesh file (.cas.h5)

    std::string dataFileName = R"(/home/s2/Documents/vm_shared/fluent4/data2/FFF.1-3-01100.dat.h5)";
    fluentReader->SetDataFileName(dataFileName.c_str());  // Set the mesh file (.cas.h5)

    fluentReader->Update(); // Update the reader

    vtkMultiBlockDataSet* set = fluentReader->GetOutput();
    spdlog::info("FluentInterpolator::LoadDataFrame; dataFile {}, blocks {}; cells {}; points {}", fileName,
                 set->GetNumberOfBlocks(), set->GetNumberOfCells(), set->GetNumberOfPoints());


    vtkUnstructuredGrid* extractedGrid = vtkUnstructuredGrid::SafeDownCast(fluentReader->GetOutput()->GetBlock(0));



    vtkCellData* cellData = extractedGrid->GetCellData();
    if (!cellData || cellData->GetNumberOfArrays() == 0) {
        std::cerr << "Error: No per-cell scalar data found!" << std::endl;
        throw std::runtime_error("celldata");
    }
    spdlog::info("celldata arrays {}", cellData->GetNumberOfArrays());

    mapper->SetScalarModeToUseCellData();  // Use cell data for coloring
//    mapper->SelectColorArray("SV_U");
    mapper->ScalarVisibilityOn();
    mapper->SetLookupTable(lut);
    lut->SetTableRange(-0.1, 0.1);

    extractedGrid->GetCellData()->SetActiveScalars("SV_U");

    //points_polydata->GetPointData()->SetActiveScalars("visualized_values");
    mapper->UseLookupTableScalarRangeOn();

//    grid->ShallowCopy(extractedGrid);


//    mapper->SetInputData(grid);
    mapper->SetInputData(extractedGrid);
    //    grid_mapper->SetLookupTable(hueLut);

    actor->SetMapper(mapper);


    actor->GetProperty()->SetEdgeVisibility(true);
    actor->GetProperty()->SetEdgeColor(0.2,0.2,0.2);
    actor->GetProperty()->LightingOff();
    actor->PickableOff();
    actor->GetProperty()->SetColor(0.1,0.1,0.1);



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

