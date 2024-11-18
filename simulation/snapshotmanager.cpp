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
    H5::DataSet dataset_grains = file.openDataSet("llGrainIDs");
    hsize_t nPoints;
    dataset_grains.getSpace().getSimpleExtentDims(&nPoints, NULL);
    model->prms.nPtsTotal = nPoints;
    spdlog::info("nPoints {}; grid [{} x {}]", nPoints,model->prms.GridXTotal, model->prms.GridY);

    // allocate space host-side
    model->gpu.hssoa.Allocate(nPoints, grid_total);
    model->gpu.hssoa.size = nPoints;

    std::vector<float> buffer(nPoints);

    // read
    file.openDataSet("x").read(buffer.data(), H5::PredType::NATIVE_FLOAT);
    t_PointReal *hssoa_ptr = model->gpu.hssoa.getPointerToPosX();
    for(int i=0;i<nPoints;i++) hssoa_ptr[i] = (t_PointReal)buffer[i];

    file.openDataSet("y").read(buffer.data(), H5::PredType::NATIVE_FLOAT);
    hssoa_ptr = model->gpu.hssoa.getPointerToPosY();
    for(int i=0;i<nPoints;i++) hssoa_ptr[i] = (t_PointReal)buffer[i];

    std::vector<uint32_t> buffer2(nPoints);
    dataset_grains.read(buffer2.data(), H5::PredType::NATIVE_UINT32);
    hssoa_ptr = model->gpu.hssoa.getPointerToLine(SimParams::idx_utility_data);
    for(int i=0;i<nPoints;i++)
    {
        *reinterpret_cast<uint32_t*>(&hssoa_ptr[i]) = buffer2[i];
    }

    // read volume attribute
    H5::Attribute att_volume = dataset_grains.openAttribute("volume");
    float volume;
    att_volume.read(H5::PredType::NATIVE_FLOAT, &volume);
    model->prms.Volume = volume;

    int offsetIncluded = 0;
    if(H5Aexists(dataset_grains.getId(), "offsetIncluded") > 0)
    {
        dataset_grains.openAttribute("offsetIncluded").read(H5::PredType::NATIVE_INT, &offsetIncluded);
    }

    // GRID
    dataset_grid.read(model->gpu.hssoa.grid_status_buffer.data(), H5::PredType::NATIVE_UINT8);

    file.close();

    // get block dimensions
    auto boundaries = model->gpu.hssoa.getBlockDimensions();
    model->prms.xmin = boundaries.first.x();
    model->prms.ymin = boundaries.first.y();
    model->prms.xmax = boundaries.second.x();
    model->prms.ymax = boundaries.second.y();
    spdlog::info("block extent x: [{}, {}], y: [{}, {}]", model->prms.xmin, model->prms.xmax, model->prms.ymin, model->prms.ymax);


    const t_PointReal &h = model->prms.cellsize;
    const t_PointReal box_x = model->prms.GridXTotal*h;
    const t_PointReal length = model->prms.xmax - model->prms.xmin;
    const t_PointReal x_offset = (box_x - length)/2;
    const t_PointReal y_offset = 2*h;

    PointVector2r offset(x_offset, y_offset);
//    if(!offsetIncluded) model->gpu.hssoa.offsetBlock(offset);
    model->gpu.hssoa.convertToIntegerCellFormat(model->prms.cellsize);
    model->gpu.hssoa.RemoveDisabledAndSort(model->prms.GridY);
    model->gpu.hssoa.InitializeBlock();

    // set indenter starting position
    const float block_left = x_offset;
    const float block_top = model->prms.ymax + y_offset;

    const float r = model->prms.IndDiameter/2;
    const float ht = r - model->prms.IndDepth;
    const float x_ind_offset = sqrt(r*r - ht*ht);

    model->prms.indenter_x = floor((block_left-x_ind_offset)/h)*h;
    if(model->prms.SetupType == 0)
        model->prms.indenter_y = block_top + ht;
    else if(model->prms.SetupType == 1)
        model->prms.indenter_y = ceil(block_top/h)*h;
    else
    {
        // for buoyancy testing
        model->prms.indenter_x = model->prms.indenter_y = -1;
    }


    model->prms.indenter_x_initial = model->prms.indenter_x;
    model->prms.indenter_y_initial = model->prms.indenter_y;

    // particle volume and mass
    model->prms.ParticleVolume = model->prms.Volume/nPoints;
    model->prms.ComputeHelperVariables();

    // allocate GPU partitions
    model->gpu.initialize();
    model->gpu.split_hssoa_into_partitions();
    model->gpu.allocate_arrays();
    model->gpu.transfer_to_device();

    model->Reset();
    model->Prepare();

    spdlog::info("LoadRawPoints done; nPoitns {}\n",nPoints);
}


