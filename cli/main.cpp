#include <iostream>
#include <functional>
#include <string>
#include <filesystem>
#include <atomic>
#include <thread>
#include <chrono>

#include <cxxopts.hpp>
#include "model.h"


int main(int argc, char** argv)
{
    // parse options
    std::string parameter_filename;
    std::string resume_filename; // Defaults to empty

    cxxopts::Options options("Ice MPM", "CLI version of MPM simulation");
    options.add_options()
        // Define options and link them to variables
        ("parameter", "Parameter file (JSON, required)", cxxopts::value<std::string>(parameter_filename))
        ("r,resume", "Resume from snapshot file <filename>", cxxopts::value<std::string>(resume_filename));

    // Mark 'parameter' as the positional argument
    options.parse_positional({"parameter"});
    auto result = options.parse(argc, argv); // Parse arguments

    // Check if required parameter file was provided
    if (!result.count("parameter")) {
        LOGV("Error: Parameter file argument is required."); // Minimal error output
        throw std::runtime_error("Parameter file is required."); // Still need to signal failure
    }

    icy::Model model;
    std::atomic<bool> request_terminate = false;


    model.LoadParameterFile(parameter_filename, resume_filename);
    model.prms.Printout();
    model.Prepare();


    model.gpu.transfer_completion_callback = [&](){
        LOGR("cycle callback {}; ", model.prms.AnimationFrameNumber());
        model.UnlockCycleMutex();
    };


    bool step_result;
    do
    {
        step_result = model.Step();
    } while(step_result);


    model.gpu.synchronize();
    std::cout << "cm done\n";
    return 0;
}
