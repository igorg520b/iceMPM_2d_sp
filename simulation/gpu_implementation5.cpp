#include "gpu_implementation5.h"
#include "parameters_sim.h"
#include "point.h"
#include "model.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include <Eigen/Core>
#include <Eigen/LU>

#include <spdlog/spdlog.h>

using namespace Eigen;


void GPU_Implementation5::reset_grid()
{
    cudaError_t err;
    for(GPU_Partition &p : partitions)
    {
        err = cudaSetDevice(p.Device);
        if(err != cudaSuccess) throw std::runtime_error("reset_grid set device");
        err = cudaEventRecord(p.event_10_cycle_start, p.streamCompute);
        if(err != cudaSuccess) throw std::runtime_error("reset_grid event");
        p.reset_grid();
    }
}

void GPU_Implementation5::p2g()
{
    cudaError_t err;

    partitions.front().p2g();
    cudaEventRecord(partitions.front().event_20_grid_halo_sent, partitions.front().streamCompute);
}



void GPU_Implementation5::update_nodes(float simulation_time, float windSpeed, float windAngle)
{
//    float interpolation_coeff = model->wind_interpolator.interpolationCoeffFromTime(simulation_time);
    float interpolation_coeff_w = model->fluent_interpolatror.position;

//    float windSpeed = std::min(0+simulation_time*(40./5.e5), 40.0);
    float alphaRad = (windAngle + 180.0) * M_PI / 180.0;

    // Compute the x and y components
    float vx = windSpeed * std::sin(alphaRad); // Eastward component
    float vy = windSpeed * std::cos(alphaRad); // Northward component

    GridVector2r vWind(vx,vy);

    for(GPU_Partition &p : partitions)
    {
        p.update_nodes(simulation_time, vWind, interpolation_coeff_w);
        cudaError_t err = cudaEventRecord(p.event_40_grid_updated, p.streamCompute);
        if(err != cudaSuccess) throw std::runtime_error("update_nodes cudaEventRecord");
    }
}

void GPU_Implementation5::g2p(const bool recordPQ, const bool enablePointTransfer, int applyGlensLaw)
{
    for(GPU_Partition &p : partitions)
    {
        p.g2p(recordPQ, enablePointTransfer, applyGlensLaw);
        cudaError_t err = cudaEventRecord(p.event_50_g2p_completed, p.streamCompute);
        if(err != cudaSuccess) throw std::runtime_error("g2p cudaEventRecord");
    }
}




void GPU_Implementation5::record_timings(const bool enablePointTransfer)
{
    for(GPU_Partition &p : partitions) p.record_timings(enablePointTransfer);
}



// ==========================================================================


void GPU_Implementation5::initialize()
{
    const int nPartitions = 1;

    // count available GPUs
    int deviceCount = 0;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);
    if (err != cudaSuccess) throw std::runtime_error("cudaGetDeviceCount error");
    if(deviceCount == 0) throw std::runtime_error("No avaialble CUDA devices");

    partitions.clear();
    partitions.resize(nPartitions);

    for(int i=0;i<nPartitions;i++)
    {
        GPU_Partition &p = partitions[i];
        p.initialize(i%deviceCount, i);
    }
}

void GPU_Implementation5::split_hssoa_into_partitions()
{
    spdlog::info("split_hssoa_into_partitions() start");

    GPU_Partition &p = partitions.front();
    p.nPts_partition = hssoa.size;
    p.GridX_partition = model->prms.GridXTotal;
    p.disabled_points_count = 0;
    spdlog::info("split_hssoa_into_partitions() done");
}


void GPU_Implementation5::allocate_host_arrays()
{
    // allocate space for points on hssoa
    model->gpu.hssoa.Allocate(model->prms.nPtsInitial);

    // points' colors (never copied to device)
    point_colors_rgb.resize(model->prms.nPtsInitial);

    // grid sizes
    const int modeled_grid_total = model->prms.GridXTotal * model->prms.GridYTotal;
    const int initial_image_total = model->prms.InitializationImageSizeX * model->prms.InitializationImageSizeY;

    // allocate grid arrays
    grid_status_buffer.resize(modeled_grid_total);
    original_image_colors_rgb.resize(3*initial_image_total);

    spdlog::info("GPU_Implementation5::allocate_host_arrays() completed");
}



void GPU_Implementation5::allocate_arrays()
{
    spdlog::info("GPU_Implementation5::allocate_arrays()");
    int GridX_size = partitions.front().GridX_partition;
    partitions.front().allocate(model->prms.nPtsInitial, GridX_size);
    spdlog::info("GPU_Implementation5::allocate_arrays() done");
}



void GPU_Implementation5::transfer_to_device()
{
    spdlog::info("GPU_Implementation: transfer_to_device()");
    int points_uploaded = 0;
    GPU_Partition &p = partitions.front();
    p.transfer_points_from_soa_to_device(hssoa, points_uploaded);
    points_uploaded += p.nPts_partition;
    spdlog::info("transfer_to_device() done; transferred points {}", points_uploaded);

    p.transfer_grid_data_to_device(this);
}



void GPU_Implementation5::transfer_from_device()
{
    unsigned offset_pts = 0;
    for(int i=0;i<partitions.size();i++)
    {
        GPU_Partition &p = partitions[i];
        int capacity_required = offset_pts + p.nPts_partition;
        if(capacity_required > hssoa.capacity)
        {
            spdlog::error("transfer_from_device(): capacity {} exceeded ({}) when transferring P {}",
                             hssoa.capacity, capacity_required, p.PartitionID);
            throw std::runtime_error("transfer_from_device capacity exceeded");
        }

        p.transfer_from_device(hssoa, offset_pts);
        offset_pts += p.nPts_partition;
    }
    hssoa.size = offset_pts;

    // wait until everything is copied to host
    for(int i=0;i<partitions.size();i++)
    {
        GPU_Partition &p = partitions[i];
        cudaSetDevice(p.Device);
        cudaStreamSynchronize(p.streamCompute);
        if(p.error_code)
        {
            spdlog::critical("P {}; error code {}", p.PartitionID, p.error_code);
            throw std::runtime_error("error code");
        }
    }

    if(transfer_completion_callback) transfer_completion_callback();
}


void GPU_Implementation5::synchronize()
{
    for(GPU_Partition &p : partitions)
    {
        cudaSetDevice(p.Device);
        cudaDeviceSynchronize();
    }
}

void GPU_Implementation5::update_constants()
{
    for(GPU_Partition &p : partitions) p.update_constants();
}

void GPU_Implementation5::reset_timings()
{
    for(GPU_Partition &p : partitions)
    {
        p.reset_timings();
    }
}

void GPU_Implementation5::update_wind_velocity_grid()
{
    spdlog::info("GPU_Implementation5::update_wind_velocity_grid()");
    for(GPU_Partition &p : partitions) p.update_wind_velocity_grid(model->wind_interpolator.grid);
}




/*

// ========================================= initialization and kernel execution

void CUDART_CB GPU_Implementation5::callback_from_stream(cudaStream_t stream, cudaError_t status, void *userData)
{
    // simulation data was copied to host memory -> proceed with processing of this data
    GPU_Implementation5 *gpu = reinterpret_cast<GPU_Implementation5*>(userData);
    // any additional processing here
    if(gpu->transfer_completion_callback) gpu->transfer_completion_callback();
}
*/
