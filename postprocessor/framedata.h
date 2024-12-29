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
    bool dataLoaded = false;
    SimParams prms;

    FrameData();
    void LoadHDF5Frame(std::string frameFileName);


    std::vector<uint8_t> grid_status_buffer;
    std::vector<uint8_t> count;
    std::vector<float> vis_r, vis_g, vis_b, vis_alpha, vis_Jpinv, vis_P, vis_Q, vis_vx, vis_vy;

private:
    void LoadWindData(std::string initialGridAndWindFileName);
};

#endif // FRAMEDATA_H
