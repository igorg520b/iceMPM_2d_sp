#ifndef GPU_PARTITION_H
#define GPU_PARTITION_H

#include <Eigen/Core>
#include <Eigen/LU>
#include <spdlog/spdlog.h>

#include <cuda.h>
#include <cuda_runtime.h>

#include <functional>
#include <vector>

#include "parameters_sim.h"
#include "point.h"
#include "host_side_soa.h"


// kernels
__global__ void partition_kernel_p2g(const int gridX, const int gridX_offset, const int pitch_grid,
                              const int count_pts, const int pitch_pts,
                                     const float *buffer_pts, float *buffer_grid);

__global__ void partition_kernel_update_nodes(const Eigen::Vector2f indCenter,
                                              const int nNodes, const int gridX_offset, const int pitch_grid,
                                              float *_buffer_grid, float *indenter_force_accumulator,
                                              float simulation_time);

__global__ void partition_kernel_g2p(const bool recordPQ, const bool enablePointTransfer,
                                     const int gridX, const int gridX_offset, const int pitch_grid,
                                     const int count_pts, const int pitch_pts,
                                     float *buffer_pts, const float *buffer_grid);




__device__ void svd(const float a[4], float u[4], float sigma[2], float v[4]);
__device__ void svd2x2(const Eigen::Matrix2f &mA, Eigen::Matrix2f &mU, Eigen::Vector2f &mS, Eigen::Matrix2f &mV);

__device__ void Wolper_Drucker_Prager(icy::Point &p);
__device__ void CheckIfPointIsInsideFailureSurface(icy::Point &p);
__device__ Eigen::Matrix2f KirchhoffStress_Wolper(const Eigen::Matrix2f &F);

__device__ void ComputeSVD(icy::Point &p, const float &kappa, const float &mu);
__device__ void ComputePQ(icy::Point &p, const float &kappa, const float &mu);
__device__ void GetParametersForGrain(short grain, float &pmin, float &pmax, float &qmax, float &beta, float &mSq, float &pmin2);

__device__ Eigen::Vector2f dev_d(Eigen::Vector2f Adiag);
__device__ Eigen::Matrix2f dev(Eigen::Matrix2f A);


struct GPU_Partition
{
    GPU_Partition();
    ~GPU_Partition();

    // preparation
    void initialize(int device, int partition);
    void allocate(int n_points_capacity, int grid_x_capacity);
    void transfer_points_from_soa_to_device(HostSideSOA &hssoa, int point_idx_offset);
    void update_constants();
    void transfer_from_device(HostSideSOA &hssoa, int point_idx_offset);

    // calculation
    void reset_grid();
    void reset_indenter_force_accumulator();
    void p2g();
    void update_nodes(float simulation_time);
    void g2p(const bool recordPQ, const bool enablePointTransfer);

    // analysis
    void reset_timings();
    void record_timings(const bool enablePointTransfer);
    void normalize_timings(int cycles);

    // host-side data
    int PartitionID;
    int Device;
    static icy::SimParams *prms;

    size_t nPtsPitch, nGridPitch; // in number of elements(!), for coalesced access on the device
    int nPts_partition;    // actual number of points (including disabled)
    int nPts_disabled;      // count the disabled points in this partition
    int GridX_partition;   // size of the portion of the grid for which this partition is "responsible"
    int GridX_offset;      // index where the grid starts in this partition

    float *host_side_indenter_force_accumulator;

    // stream and events
    cudaStream_t streamCompute;

    cudaEvent_t event_10_cycle_start;
    cudaEvent_t event_20_grid_halo_sent;
    cudaEvent_t event_30_halo_accepted;
    cudaEvent_t event_40_grid_updated;
    cudaEvent_t event_50_g2p_completed;
    cudaEvent_t event_70_pts_sent;
    cudaEvent_t event_80_pts_accepted;

    bool initialized = false;
    uint8_t error_code = 0;

    // device-side data
    float *pts_array, *grid_array, *indenter_force_accumulator;

    // frame analysis
    float timing_10_P2GAndHalo;
    float timing_20_acceptHalo;
    float timing_30_updateGrid;
    float timing_40_G2P;
    float timing_60_ptsSent;
    float timing_70_ptsAccepted;
    float timing_stepTotal;
};


#endif // GPU_PARTITION_H
