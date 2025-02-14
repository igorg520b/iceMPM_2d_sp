#include "model.h"
#include <spdlog/spdlog.h>
#include <algorithm>
#include <thread>


bool icy::Model::Step()
{
//    if(prms.SimulationStep == 0 && prms.SaveSnapshots) SaveFrameRequest(prms.SimulationStep, prms.SimulationTime);

    std::cout << '\n';
    spdlog::info("step {} ({}) started; sim_time {:>6.3}; host pts {}; cap {}",
                 prms.SimulationStep, prms.AnimationFrameNumber(), prms.SimulationTime, gpu.hssoa.size, gpu.hssoa.capacity);

    int count_unupdated_steps = 0;
    gpu.reset_timings();
    const float pp = 1000;
    const double st = prms.SimulationTime;
    std::pair<float, float> spd_and_angle = {
                                             10+(20+st/pp)*sin(SimParams::pi * st/pp),
                                             370+10*sin(SimParams::pi * st/(10*60.))};
    double simulation_time;
    do
    {
        const int step = prms.SimulationStep + count_unupdated_steps;
        simulation_time = prms.InitialTimeStep * step;

        gpu.reset_grid();
        gpu.p2g();

//        if(prms.UseWindData && wind_interpolator.setTime(simulation_time)) gpu.update_wind_velocity_grid();
//        if(prms.UseCurrentData && fluent_interpolatror.SetTime(simulation_time))
//            gpu.partitions.front().update_water_flow_grid(fluent_interpolatror.getFramePtr(0,0),
 //                                                         fluent_interpolatror.getFramePtr(0,1),
  //                                                        fluent_interpolatror.getFramePtr(1,0),
  //                                                        fluent_interpolatror.getFramePtr(1,1));

        if(prms.UseCurrentData && fluent_interpolatror.SetTime(0))
            gpu.partitions.front().update_water_flow_grid(fluent_interpolatror.getFramePtr(0,0),
                                                          fluent_interpolatror.getFramePtr(0,1),
                                                          fluent_interpolatror.getFramePtr(1,0),
                                                          fluent_interpolatror.getFramePtr(1,1));


        gpu.update_nodes(simulation_time, spd_and_angle.first, spd_and_angle.second);
        const bool isCycleStart = step % prms.UpdateEveryNthStep == 0;
        gpu.g2p(isCycleStart, false, false);
        gpu.record_timings(false);

        count_unupdated_steps++;
        if(intentionalSlowdown)
        {
            gpu.synchronize();
            std::this_thread::sleep_for(std::chrono::milliseconds(intentionalSlowdown));
        }
    } while((prms.SimulationStep+count_unupdated_steps) % prms.UpdateEveryNthStep != 0);

    processing_current_cycle_data.lock();   // if locked, previous results are not yet processed by the host
    accessing_point_data.lock();

    gpu.transfer_from_device();
    spdlog::info("finished {} ({}); host pts {}; cap {}; windSpeed {}; windBearing {}", prms.SimulationEndTime,
                 prms.AnimationFrameNumber(), gpu.hssoa.size, gpu.hssoa.capacity, spd_and_angle.first, spd_and_angle.second);
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
        gpu.hssoa.RemoveDisabledAndSort(prms.GridYTotal);
        gpu.split_hssoa_into_partitions();
        gpu.transfer_to_device();
        SyncTopologyRequired = true;
        spdlog::info("Model::Step() rebalancing done");
    }

    accessing_point_data.unlock();
    if(prms.SaveSnapshots) SaveFrameRequest(prms.SimulationStep, prms.SimulationTime);


    return (prms.SimulationTime < prms.SimulationEndTime);
}




icy::Model::Model() : frame_ready(false), done(false)
{
    fluent_interpolatror.gpu = &gpu;
    fluent_interpolatror.prms = &prms;
    snapshot.model = this;
    prms.SimulationStep = 0;
    prms.SimulationTime = 0;
    SyncTopologyRequired = true;

    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("logs/multisink.txt", true);
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    auto lg = std::make_shared<spdlog::logger>("multi_sink", spdlog::sinks_init_list({console_sink, file_sink}));
    spdlog::set_default_logger(lg);
    spdlog::set_pattern("%v");

    prms.Reset();
    gpu.model = this;
    GPU_Partition::prms = &this->prms;
    saver_thread = std::thread(&icy::Model::SaveThread, this);
    spdlog::info("Model constructor");
}


void icy::Model::SaveFrameRequest(int SimulationStep, double SimulationTime)
{
    {
        std::lock_guard<std::mutex> lock(frame_mutex);
        spdlog::info("icy::Model::SaveFrameRequest; step {}",SimulationStep);
        saving_SimulationStep = SimulationStep;
        saving_SimulationTime = SimulationTime;
        frame_ready = true; // Indicate that a new frame is ready
    }
    frame_cv.notify_one(); // Notify the saver thread
}


icy::Model::~Model()
{
    {
    std::lock_guard<std::mutex> lock(frame_mutex);
    done = true; // Signal that we're done
    saving_SimulationStep = -1;
    }
    frame_cv.notify_one(); // Notify the saver thread
    saver_thread.join();   // Wait for the thread to finish
}


// Frame-saving thread function
void icy::Model::SaveThread() {
    while (true) {
        // Wait for a frame to save or for the simulation to finish
        {
            std::unique_lock<std::mutex> lock(frame_mutex);
            frame_cv.wait(lock, [this] { return frame_ready || done.load(); });

            if (frame_ready) {
                frame_ready = false; // Mark the frame as consumed
            } else if (done.load()) {
                break; // Exit if simulation is finished
            }
        }

        // Save the frame outside the critical section
        if (saving_SimulationStep != -1) {
            bool saveSnapshot = (saving_SimulationStep/prms.UpdateEveryNthStep)%SimParams::snapshotFrequency == 0;
            accessing_point_data.lock();
            if(saveSnapshot) snapshot.SaveSnapshot(saving_SimulationStep, saving_SimulationTime);
            snapshot.SaveFrame(saving_SimulationStep, saving_SimulationTime);
            accessing_point_data.unlock();

            saving_SimulationStep = -1;
        }
    }
}



void icy::Model::UnlockCycleMutex()
{
    // current data was handled by host - allow next cycle to proceed
    processing_current_cycle_data.unlock();
}


void icy::Model::Prepare()
{
    spdlog::info("icy::Model::Prepare()");
    //abortRequested = false;
    gpu.update_constants();

    fluent_interpolatror.SetTime(prms.SimulationTime);
    gpu.partitions.front().update_water_flow_grid(fluent_interpolatror.getFramePtr(0,0),
                                                  fluent_interpolatror.getFramePtr(0,1),
                                                  fluent_interpolatror.getFramePtr(1,0),
                                                  fluent_interpolatror.getFramePtr(1,1));
}
