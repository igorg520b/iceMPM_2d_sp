#ifndef FRAMEDATA_H
#define FRAMEDATA_H

#include <string>
#include <vector>
#include <algorithm>
#include <mutex>
#include <iostream>

#include <H5Cpp.h>
#include <spdlog/spdlog.h>
#include <Eigen/Core>
#include "parameters_sim.h"
#include "vtk_visualization.h"
#include "generalgriddata.h"

#include <vtkRenderWindow.h>
#include <vtkWindowToImageFilter.h>
#include <vtkJPEGWriter.h>
#include <vtkRenderer.h>
#include <vtkOpenGLRenderWindow.h>
#include <vtkCamera.h>



struct FrameData
{
    FrameData(GeneralGridData &ggd_);

    GeneralGridData &ggd;
    VTKVisualization representation;
    double SimulationTime;

    std::vector<uint8_t> count;
    std::vector<float> vis_Jpinv, vis_P, vis_Q, vis_vx, vis_vy;
    std::vector<uint8_t> rgb;

    void LoadHDF5Frame(int frameNumber);

    void SetUpOffscreenRender(const FrameData &guiFD, vtkCamera* sourceGuiCamera);
    void RenderFrame(VTKVisualization::VisOpt visopt);

    vtkNew<vtkRenderer> renderer;
    int frame = 0;
private:
    int SimulationStep;


    // offscreenrender
    vtkNew<vtkRenderWindow> offscreenRenderWindow;
    vtkNew<vtkWindowToImageFilter> windowToImageFilter;
    vtkNew<vtkJPEGWriter> writer;
};

#endif // FRAMEDATA_H
