#include "model.h"
#include <spdlog/spdlog.h>
#include <algorithm>


icy::Model::Model()
{
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("logs/multisink.txt", true);
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    auto lg = std::make_shared<spdlog::logger>("multi_sink", spdlog::sinks_init_list({console_sink, file_sink}));
    spdlog::set_default_logger(lg);
    spdlog::set_pattern("%v");

    prms.Reset();
    gpu.model = this;
    GPU_Partition::prms = &this->prms;
    spdlog::info("Model constructor");
}

icy::Model::~Model() {}

bool icy::Model::Step()
{
    float simulation_time = prms.SimulationTime;
    std::cout << '\n';
    spdlog::info("step {} ({}) started; sim_time {:>6.3}; host pts {}; cap {}",
                 prms.SimulationStep, prms.AnimationFrameNumber(), simulation_time,
                 gpu.hssoa.size, gpu.hssoa.capacity);

    int count_unupdated_steps = 0;
    gpu.reset_timings();

    do
    {
        simulation_time += prms.InitialTimeStep;
        const int step = prms.SimulationStep+count_unupdated_steps;

        gpu.reset_grid();
        gpu.p2g();
        gpu.update_nodes(simulation_time);
        const bool isZeroStep = step % prms.UpdateEveryNthStep == 0;
        gpu.g2p(isZeroStep, false, 1);
        gpu.record_timings(false);

        count_unupdated_steps++;
    } while((prms.SimulationStep+count_unupdated_steps) % prms.UpdateEveryNthStep != 0);

    processing_current_cycle_data.lock();   // if locked, previous results are not yet processed by the host
    accessing_point_data.lock();

    gpu.transfer_from_device();
    spdlog::info("finished {} ({}); host pts {}; cap {}", prms.SimulationEndTime,
                 prms.AnimationFrameNumber(), gpu.hssoa.size, gpu.hssoa.capacity);
    prms.SimulationTime = simulation_time;
    prms.SimulationStep += count_unupdated_steps;

    // print out timings
    spdlog::info("{:^3s} {:^8s} {:^8s} {:^7s} | {:^5s} {:^5s} {:^5s} | {:^5s} {:^5s} {:^5s} {:^5s} {:^5s} | {:^6s}",
                 "P-D",  "pts", "free", "dis",  "p2g",  "s2",  "S12",     "u",  "g2p", "psnt", "prcv","S36", "tot");
    GPU_Partition &p = gpu.partitions.front();
    p.normalize_timings(count_unupdated_steps);
    spdlog::info("{:>1}-{:>1} {:>8} {:>8} {:>7} | {:>5.1f} {:>5.1f} {:>5.1f} | {:>5.1f} {:>5.1f} {:>5.1f} {:5.1f} {:5.1f} | {:>6.1f}",
                 p.PartitionID, p.Device, p.nPts_partition, (p.nPtsPitch-p.nPts_partition), p.disabled_points_count,
                 p.timing_10_P2GAndHalo, p.timing_20_acceptHalo, (p.timing_10_P2GAndHalo + p.timing_20_acceptHalo),
                 p.timing_30_updateGrid, p.timing_40_G2P, p.timing_60_ptsSent, p.timing_70_ptsAccepted,
                 (p.timing_30_updateGrid + p.timing_40_G2P + p.timing_60_ptsSent + p.timing_70_ptsAccepted),
                 p.timing_stepTotal);

    const float disabled_proportion = (float)p.disabled_points_count/p.nPts_partition;
    if(disabled_proportion > SimParams::disabled_pts_proportion_threshold)
    {
        spdlog::info("Model::Step() squeezing and sorting HSSOA");
        gpu.hssoa.RemoveDisabledAndSort(prms.GridY);
        gpu.split_hssoa_into_partitions();
        gpu.transfer_to_device();
        SyncTopologyRequired = true;
        spdlog::info("Model::Step() rebalancing done");
    }

    accessing_point_data.unlock();

    return (prms.SimulationTime < prms.SimulationEndTime);
}


void icy::Model::UnlockCycleMutex()
{
    // current data was handled by host - allow next cycle to proceed
    processing_current_cycle_data.unlock();
}


void icy::Model::Reset()
{
    spdlog::info("icy::Model::Reset()");

    prms.SimulationStep = 0;
    prms.SimulationTime = 0;
    SyncTopologyRequired = true;
}

void icy::Model::Prepare()
{
    spdlog::info("icy::Model::Prepare()");
    abortRequested = false;
    gpu.update_constants();
}

