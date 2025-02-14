#ifndef FRAMEDATA_H
#define FRAMEDATA_H

#include <string>
#include <vector>
#include <algorithm>
#include <mutex>

#include <H5Cpp.h>
#include <spdlog/spdlog.h>
#include <Eigen/Core>
#include "parameters_sim.h"
#include "vtk_visualization.h"

//#include "windinterpolator.h"
//#include "fluentinterpolator.h"

#include <vtkRenderWindow.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkRenderer.h>
#include <vtkOpenGLRenderWindow.h>
#include <vtkCamera.h>

struct GeneralGridData
{
    void ScanDirectory(std::string frameFileName);
    SimParams prms;
    int countFrames;
    std::string frameDirectory;

    std::vector<uint8_t> grid_status_buffer;
    std::vector<uint8_t> original_image_colors_rgb;
};


struct FrameData
{
    FrameData();
    FrameData(const FrameData &other);

    VTKVisualization representation;

    void SetUpOffscreenRender(FrameData &guiFD, double camData[10]);
    void LoadHDF5Frame(int frameNumber);
    void LoadHDF5Frame(std::string fileName);
    void RenderFrame(VTKVisualization::VisOpt visopt, std::string_view outputDirectory);

    GeneralGridData *ggd = nullptr;
    bool dataLoaded = false;
    int SimulationStep;
    double SimulationTime;
    int frame = 0;
    std::mutex *mutex;

    std::vector<uint8_t> count;
    std::vector<float> vis_Jpinv, vis_P, vis_Q, vis_vx, vis_vy;
    std::vector<uint8_t> rgb;

//    vtkNew<vtkRenderWindow> offscreenRenderWindow;
    vtkNew<vtkRenderWindow> offscreenRenderWindow;
    vtkNew<vtkWindowToImageFilter> windowToImageFilter;
    vtkNew<vtkPNGWriter> writerPNG;
    vtkNew<vtkRenderer> renderer;
};

#endif // FRAMEDATA_H
