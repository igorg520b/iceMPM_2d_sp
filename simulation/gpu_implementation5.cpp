#include "gpu_implementation5.h"
#include "parameters_sim.h"
#include "point.h"
#include "model.h"

#include <iostream>
#include <vector>
#include <algorithm>

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
    const int &gridY = model->prms.GridY;
    const int &gridXT = model->prms.GridXTotal;

    partitions.front().p2g();

    cudaEventRecord(partitions.front().event_20_grid_halo_sent, partitions.front().streamCompute);
}



void GPU_Implementation5::update_nodes(float simulation_time)
{
    for(GPU_Partition &p : partitions)
    {
        p.update_nodes(simulation_time);
        cudaError_t err = cudaEventRecord(p.event_40_grid_updated, p.streamCompute);
        if(err != cudaSuccess) throw std::runtime_error("update_nodes cudaEventRecord");
    }
}

void GPU_Implementation5::g2p(const bool recordPQ, const bool enablePointTransfer)
{
    for(GPU_Partition &p : partitions)
    {
        p.g2p(recordPQ, enablePointTransfer);
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
    const t_PointReal &hinv = model->prms.cellsize_inv;
    const int &GridXTotal = model->prms.GridXTotal;

    GPU_Partition &p = partitions.front();
    p.nPts_disabled = 0;
    p.nPts_partition = hssoa.size;
    p.GridX_offset = 0;
    p.GridX_partition = GridXTotal;
    spdlog::info("split_hssoa_into_partitions() done");
}


void GPU_Implementation5::allocate_arrays()
{
    spdlog::info("GPU_Implementation5::allocate_arrays()");
    int GridX_size = partitions.front().GridX_partition;
    partitions.front().allocate(model->prms.nPtsTotal, GridX_size);
    spdlog::info("GPU_Implementation5::allocate_arrays() done");
}



void GPU_Implementation5::transfer_ponts_to_device()
{
    spdlog::info("GPU_Implementation: transfer_to_device()");
    int points_uploaded = 0;
    for(GPU_Partition &p : partitions)
    {
        p.transfer_points_from_soa_to_device(hssoa, points_uploaded);
        points_uploaded += p.nPts_partition;
    }
    spdlog::info("transfer_ponts_to_device() done; transferred points {}", points_uploaded);
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
    indenter_force.setZero();
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

        for(int j=0; j<model->prms.n_indenter_subdivisions; j++)
        {
            indenter_force.x() += p.host_side_indenter_force_accumulator[j*2+0];
            indenter_force.y() += p.host_side_indenter_force_accumulator[j*2+1];
        }
    }

    int count = 0;
    for(int i=0;i<hssoa.size;i++)
    {
        SOAIterator s = hssoa.begin()+i;
        if(s->getDisabledStatus()) continue;
        count++;
    }

    if(count != model->prms.nPtsTotal)
    {
        spdlog::error("tranfer: hssoa.size {}; nPts {}, count_active {}", hssoa.size, model->prms.nPtsTotal, count);

        unsigned offset_pts = 0;
        for(int i=0;i<partitions.size();i++)
        {
            GPU_Partition &p = partitions[i];
            int count_disabled_soa = 0;
            for(int i=offset_pts; i<offset_pts+p.nPts_partition; i++)
            {
                SOAIterator s = hssoa.begin()+i;
                if(s->getDisabledStatus())
                {
                    std::cout << i << ' ';
                    count_disabled_soa++;
                }
            }
            std::cout << '\n';
            offset_pts += p.nPts_partition;
            spdlog::error("P{}: size {}; disabled {}; disabled_soa {}",
                          p.PartitionID, p.nPts_partition, p.nPts_disabled, count_disabled_soa);
        }


        throw std::runtime_error("transfer_from_device(): active point count mismatch");
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

void GPU_Implementation5::reset_indenter_force_accumulator()
{
    for(GPU_Partition &p : partitions)
    {
        p.reset_indenter_force_accumulator();
        p.reset_timings();
    }
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
