#ifndef FRAMEDATA_H
#define FRAMEDATA_H

#include <string.h>
#include <H5Cpp.h>
#include <spdlog/spdlog.h>
#include <Eigen/Core>
#include "simulation/parameters_sim.h"


class FrameData
{
public:
    SimParams prms;

    FrameData();
    void LoadHDF5Frame(std::string frameFileName);


    std::vector<uint8_t> grid_status_buffer;
private:
    void LoadWindData(std::string initialGridAndWindFileName);
};

#endif // FRAMEDATA_H
