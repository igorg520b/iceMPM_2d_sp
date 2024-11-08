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

    H5::DataSet dataset_grains = file.openDataSet("llGrainIDs");
    hsize_t nPoints;
    dataset_grains.getSpace().getSimpleExtentDims(&nPoints, NULL);
    model->prms.nPtsTotal = nPoints;

    // allocate space host-side
    model->gpu.hssoa.Allocate(nPoints);
    model->gpu.hssoa.size = nPoints;

    // read
    file.openDataSet("x").read(model->gpu.hssoa.getPointerToPosX(), H5::PredType::NATIVE_FLOAT);
    file.openDataSet("y").read(model->gpu.hssoa.getPointerToPosY(), H5::PredType::NATIVE_FLOAT);
    dataset_grains.read(model->gpu.hssoa.host_buffer, H5::PredType::NATIVE_UINT32);

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

    file.close();

    // get block dimensions
    auto boundaries = model->gpu.hssoa.getBlockDimensions();
    model->prms.xmin = boundaries.first.x();
    model->prms.ymin = boundaries.first.y();
    model->prms.xmax = boundaries.second.x();
    model->prms.ymax = boundaries.second.y();
    spdlog::info("block extent x: [{}, {}], y: [{}, {}]", model->prms.xmin, model->prms.xmax, model->prms.ymin, model->prms.ymax);


    const float &h = model->prms.cellsize;
    const float box_x = model->prms.GridXTotal*h;
    const float length = model->prms.xmax - model->prms.xmin;
    const float x_offset = (box_x - length)/2;
    const float y_offset = 2*h;

    Eigen::Vector2f offset(x_offset, y_offset);
    if(!offsetIncluded) model->gpu.hssoa.offsetBlock(offset);
    model->gpu.hssoa.RemoveDisabledAndSort(model->prms.cellsize_inv, model->prms.GridY);
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
    model->gpu.transfer_ponts_to_device();

    model->Reset();
    model->Prepare();

    spdlog::info("LoadRawPoints done; nPoitns {}\n",nPoints);
}


