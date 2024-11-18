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
__global__ void partition_kernel_p2g(const int gridX, const int pitch_grid,
                              const int count_pts, const int pitch_pts,
                                     const t_PointReal *buffer_pts, t_GridReal *buffer_grid);

__global__ void partition_kernel_update_nodes(const GridVector2r indCenter,
                                              const int nNodes, const int pitch_grid,
                                              t_GridReal *_buffer_grid, t_PointReal *indenter_force_accumulator,
                                              t_PointReal simulation_time, const uint8_t *grid_status);

__global__ void partition_kernel_g2p(const bool recordPQ,
                                     const int pitch_grid,
                                     const int count_pts, const int pitch_pts,
                                     t_PointReal *buffer_pts, const t_GridReal *buffer_grid);




__device__ void svd2x2(const PointMatrix2r &mA, PointMatrix2r &mU, PointVector2r &mS, PointMatrix2r &mV);

__device__ void Wolper_Drucker_Prager(const t_PointReal &p_tr, const t_PointReal &q_tr, const t_PointReal &Je_tr,
                                      const PointMatrix2r &U, const PointMatrix2r &V, const PointVector2r &vSigmaSquared, const PointVector2r &v_s_hat_tr,
                                      PointMatrix2r &Fe, t_PointReal &Jp_inv);

__device__ void CheckIfPointIsInsideFailureSurface(uint32_t &utility_data, const uint16_t &grain,
                                                   const t_PointReal &p, const t_PointReal &q);

__device__ void ComputeSVD(const PointMatrix2r &Fe, PointMatrix2r &U, PointVector2r &vSigma, PointMatrix2r &V,
                           PointVector2r &vSigmaSquared, PointVector2r &v_s_hat_tr,
                           const t_PointReal &kappa, const t_PointReal &mu, const t_PointReal &Je_tr);

__device__ void ComputePQ(t_PointReal &Je_tr, t_PointReal &p_tr, t_PointReal &q_tr,
                          const t_PointReal &kappa, const t_PointReal &mu, const PointMatrix2r &F);

__device__ void GetParametersForGrain(int grain, t_PointReal &pmin, t_PointReal &pmax, t_PointReal &qmax,
                                      t_PointReal &beta, t_PointReal &mSq, t_PointReal &pmin2);

__device__ PointMatrix2r KirchhoffStress_Wolper(const PointMatrix2r &F);

__device__ PointVector2r dev_d(PointVector2r Adiag);

__device__ PointMatrix2r dev(PointMatrix2r A);

__device__ void CalculateWeightCoeffs(const PointVector2r &pos, PointArray2r ww[3]);

__device__ PointMatrix2r Water(const t_PointReal J);


struct GPU_Partition
{
    GPU_Partition();
    ~GPU_Partition();

    // preparation
    void initialize(int device, int partition);
    void allocate(int n_points_capacity, int grid_x_capacity);
    void transfer_points_from_soa_to_device(HostSideSOA &hssoa, int point_idx_offset);
    void transfer_grid_data_to_device(HostSideSOA &hssoa);
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
    static SimParams *prms;

    size_t nPtsPitch, nGridPitch; // in number of elements(!), for coalesced access on the device
    int nPts_partition;    // actual number of points (including disabled)
    int nPts_disabled;      // count the disabled points in this partition
    int GridX_partition;   // size of the portion of the grid for which this partition is "responsible"
    int GridX_offset;      // index where the grid starts in this partition

    t_PointReal *host_side_indenter_force_accumulator;

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
    t_PointReal *pts_array, *indenter_force_accumulator;
    t_GridReal *grid_array;
    uint8_t *grid_status_array;

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
