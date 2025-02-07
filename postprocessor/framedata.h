#ifndef FRAMEDATA_H
#define FRAMEDATA_H

#include <string>
#include <vector>

#include <H5Cpp.h>
#include <spdlog/spdlog.h>
#include <Eigen/Core>
#include "parameters_sim.h"
//#include "windinterpolator.h"
//#include "fluentinterpolator.h"

class FrameData
{
public:
    bool dataLoaded = false;
    SimParams prms;
//    WindInterpolator windInterpolator;
//    FluentInterpolator fluentInterpolator;
    std::string frameDirectory;
    std::vector<bool> availableFrames;

    FrameData();
    void LoadHDF5Frame(std::string frameFileName, bool loadGridAndWind = true);
    void ScanDirectory(std::string frameFileName);

    std::vector<uint8_t> grid_status_buffer;
    std::vector<uint8_t> original_image_colors_rgb;

    std::vector<uint8_t> count;
    std::vector<float> vis_Jpinv, vis_P, vis_Q, vis_vx, vis_vy;
    std::vector<uint8_t> rgb;

private:
};

#endif // FRAMEDATA_H
