#include <filesystem>
#include <regex>
#include "framedata.h"

FrameData::FrameData() {}


void FrameData::ScanDirectory(std::string frameFileName)
{
    std::filesystem::path framePath(frameFileName);

    // Navigate to the target file by appending the relative path
    std::filesystem::path gridAndWindPath = framePath.parent_path();
    frameDirectory = gridAndWindPath.string();


    // Read the contents of the directory
    for (const auto& entry : std::filesystem::directory_iterator(frameDirectory)) {
        if (entry.is_regular_file()) {
            std::string filename = entry.path().filename().string();

            // Check if the filename matches the pattern "frame_XXXXX.h5"
            std::regex pattern(R"(frame_(\d{5})\.h5)");
            std::smatch match;
            if (std::regex_match(filename, match, pattern)) {
                // Extract the five-digit integer from the filename
                int frameIndex = std::stoi(match[1].str());

                // Ensure the availableFrames vector is large enough
                if (frameIndex >= availableFrames.size()) {
                    availableFrames.resize(frameIndex + 1, false);
                }

                // Set the corresponding element to true
                availableFrames[frameIndex] = true;
            }
        }
    }

    spdlog::info("availableFrames size {}",availableFrames.size());
    spdlog::info("frameDirectory {}", frameDirectory);
}



void FrameData::LoadHDF5Frame(std::string frameFileName, bool loadGridAndWind)
{
    spdlog::info("LoadHDF5Frame: {}; {}",frameFileName, loadGridAndWind);
    if(loadGridAndWind)
    {
        // Convert to std::filesystem::path for manipulation
        std::filesystem::path framePath(frameFileName);

        // Navigate to the target file by appending the relative path
        std::filesystem::path gridAndWindPath = framePath.parent_path() / "../snapshots/_grid_and_wind.h5";

        // Normalize the path (resolving ".." and ".")
        gridAndWindPath = std::filesystem::canonical(gridAndWindPath);

        // Convert back to std::string if needed
        std::string resultPath = gridAndWindPath.string();

        // open the initial data file
        spdlog::info("LoadHDF5Frame resultPath = {}", resultPath);
        H5::H5File file(resultPath, H5F_ACC_RDONLY);

        H5::DataSet gridDataset = file.openDataSet("grid");
        prms.ReadParametersFromHDF5Attributes(gridDataset);
        prms.Printout();

        // Get the dataspace of the dataset
        H5::DataSpace dataspace = gridDataset.getSpace();

        // Get the dimensions of the dataset
        hsize_t dims[2];
        dataspace.getSimpleExtentDims(dims, nullptr);
        size_t gridXTotal = dims[0];
        size_t gridY = dims[1];

        if(gridXTotal != prms.GridXTotal || gridY != prms.GridY)
            throw std::runtime_error("LoadHDF5Frame grid size mismatch");

        // Resize the vector to hold the data
        grid_status_buffer.resize(gridXTotal * gridY);

        // Read the dataset into the vector
        gridDataset.read(grid_status_buffer.data(), H5::PredType::NATIVE_UINT8);

        windInterpolator.ReadFromOwnHDF5(file);
    }

    size_t gridSize = prms.GridXTotal*prms.GridY;
    count.resize(gridSize);
    vis_r.resize(gridSize);
    vis_g.resize(gridSize);
    vis_b.resize(gridSize);
    vis_alpha.resize(gridSize);
    vis_vx.resize(gridSize);
    vis_vy.resize(gridSize);

    vis_Jpinv.resize(gridSize);
    vis_P.resize(gridSize);
    vis_Q.resize(gridSize);


    H5::H5File file(frameFileName, H5F_ACC_RDONLY);
    H5::DataSet ds_count = file.openDataSet("count");

    ds_count.openAttribute("SimulationStep").read(H5::PredType::NATIVE_INT, &prms.SimulationStep);
    ds_count.openAttribute("SimulationTime").read(H5::PredType::NATIVE_DOUBLE, &prms.SimulationTime);

    hsize_t dims[2];
    ds_count.getSpace().getSimpleExtentDims(dims, nullptr);
    if(dims[0] != prms.GridXTotal || dims[1] != prms.GridY) throw std::runtime_error("LoadHDF5Frame frame grid size mismatch");
    ds_count.read(count.data(), H5::PredType::NATIVE_UINT8);

    file.openDataSet("vis_vx").read(vis_vx.data(), H5::PredType::NATIVE_FLOAT);
    file.openDataSet("vis_vy").read(vis_vy.data(), H5::PredType::NATIVE_FLOAT);
    file.openDataSet("vis_r").read(vis_r.data(), H5::PredType::NATIVE_FLOAT);
    file.openDataSet("vis_g").read(vis_g.data(), H5::PredType::NATIVE_FLOAT);
    file.openDataSet("vis_b").read(vis_b.data(), H5::PredType::NATIVE_FLOAT);
    file.openDataSet("vis_alpha").read(vis_alpha.data(), H5::PredType::NATIVE_FLOAT);

    file.openDataSet("vis_Jpinv").read(vis_Jpinv.data(), H5::PredType::NATIVE_FLOAT);
    file.openDataSet("vis_P").read(vis_P.data(), H5::PredType::NATIVE_FLOAT);
    file.openDataSet("vis_Q").read(vis_Q.data(), H5::PredType::NATIVE_FLOAT);

    dataLoaded = true;
}


void FrameData::LoadWindData(std::string initialGridAndWindFileName)
{
    H5::H5File file(initialGridAndWindFileName, H5F_ACC_RDONLY);

}
