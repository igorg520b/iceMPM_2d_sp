#include <filesystem>
#include "framedata.h"

FrameData::FrameData() {}


void FrameData::LoadHDF5Frame(std::string frameFileName)
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


}


void FrameData::LoadWindData(std::string initialGridAndWindFileName)
{
    H5::H5File file(initialGridAndWindFileName, H5F_ACC_RDONLY);

}
