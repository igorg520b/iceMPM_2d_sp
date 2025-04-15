#include <iostream>
#include <functional>
#include <string>
#include <filesystem>
#include <atomic>
#include <thread>
#include <chrono>

#include <cxxopts.hpp>
#include <spdlog/spdlog.h>

#include "model.h"
#include "snapshotmanager.h"
#include "gpu_partition.h"



int main(int argc, char** argv)
{
    icy::Model model;
    std::atomic<bool> request_terminate = false;


    // initialize the model


    // parse options
    cxxopts::Options options("Ice MPM", "CLI version of MPM simulation");

    options.add_options()
        ("file", "Snapshot file", cxxopts::value<std::string>())
        ;
    options.parse_positional({"file"});

    auto option_parse_result = options.parse(argc, argv);

    if(option_parse_result.count("file"))
    {
        std::string fileName = option_parse_result["file"].as<std::string>();
        model.snapshot.ReadSnapshot(fileName);
        model.Prepare();
    }
    else
    {
        LOGV("snapshot file must be provided");
        throw std::runtime_error("no snapshot file");
    }


    model.gpu.transfer_completion_callback = [&](){
        LOGR("cycle callback {}; ", model.prms.AnimationFrameNumber());
        model.UnlockCycleMutex();
    };
    // ensure that the folder exists

    // start the simulation thread
    std::thread simulation_thread([&](){
        bool result;
        do
        {
            result = model.Step();
        } while(!request_terminate && result);
    });

    // accept console input
    do
    {
        std::string user_input;
        std::cin >> user_input;

        if(user_input[0]=='q'){
            request_terminate = true;
            LOGV("requested to save the snapshot and terminate");
        }
    } while(!request_terminate);

    model.gpu.synchronize();

    std::cout << "cm done\n";

    return 0;
}
