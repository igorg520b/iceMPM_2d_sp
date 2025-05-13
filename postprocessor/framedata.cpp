#include "framedata.h"
#include <fmt/format.h>
#include <fmt/std.h>

#include <vtkVersionMacros.h>
#include <QDebug>

namespace fs = std::filesystem;


FrameData::FrameData(GeneralGridData &ggd_) : ggd(ggd_), representation(*this)
{
    std::cout << "VTK Version (from macros):" << std::endl;
    std::cout << "  Full Version String (VTK_VERSION): " << VTK_VERSION << std::endl;
    std::cout << "  Major Version (VTK_MAJOR_VERSION): " << VTK_MAJOR_VERSION << std::endl;
    std::cout << "  Minor Version (VTK_MINOR_VERSION): " << VTK_MINOR_VERSION << std::endl;
    std::cout << "  Build Version (VTK_BUILD_VERSION): " << VTK_BUILD_VERSION << std::endl;
}

void FrameData::SetUpOffscreenRender(const FrameData &guiFD, vtkCamera* sourceGuiCamera)
{
    std::copy(std::begin(guiFD.representation.ranges), std::end(guiFD.representation.ranges), std::begin(this->representation.ranges));

    offscreenRenderWindow->Initialize();
    offscreenRenderWindow->SetOffScreenRendering(true);
    offscreenRenderWindow->SetSize(1920, 1080);
    offscreenRenderWindow->DoubleBufferOff();

    // renderer -> offscreen
    offscreenRenderWindow->AddRenderer(renderer);
    renderer->SetBackground(1.0, 1.0, 1.0);

    renderer->AddActor(representation.raster_actor);
    renderer->AddActor(representation.actor_text);
    renderer->AddActor(representation.scalarBar);
    renderer->AddActor(representation.actor_text_title);

    windowToImageFilter->SetInput(offscreenRenderWindow);
    windowToImageFilter->SetScale(1);
    windowToImageFilter->SetInputBufferTypeToRGB();
//    windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
    windowToImageFilter->ReadFrontBufferOn();
    writer->SetInputConnection(windowToImageFilter->GetOutputPort());


    // set up camera
    renderer->ResetCamera();
    vtkCamera* cam = renderer->GetActiveCamera();
    cam->DeepCopy(sourceGuiCamera);
    cam->ParallelProjectionOn();
    cam->Modified();
}


void FrameData::LoadHDF5Frame(int frameNumber)
{
    this->frame = frameNumber;
    //f00001.h5
    std::string fileName = fmt::format(fmt::runtime("{}/f{:05d}.h5"), ggd.frameDirectory, frame);

    snapshot.LoadFrameCompressed(fileName, SimulationStep, SimulationTime);


    /*
    size_t gridSize = ggd.prms.GridXTotal*ggd.prms.GridYTotal;
    count.resize(gridSize);
    vis_vx.resize(gridSize);
    vis_vy.resize(gridSize);
    vis_Jpinv.resize(gridSize);
    vis_P.resize(gridSize);
    vis_Q.resize(gridSize);
    rgb.resize(gridSize*3);

    LOGR("FrameData::LoadHDF5Frame {}",fileName);
    H5::H5File file(fileName, H5F_ACC_RDONLY);

    H5::DataSet ds_count = file.openDataSet("count");
    ds_count.openAttribute("SimulationStep").read(H5::PredType::NATIVE_INT, &SimulationStep);
    ds_count.openAttribute("SimulationTime").read(H5::PredType::NATIVE_DOUBLE, &SimulationTime);

    hsize_t dims[2];
    ds_count.getSpace().getSimpleExtentDims(dims, nullptr);
    if(dims[0] != ggd.prms.GridXTotal || dims[1] != ggd.prms.GridYTotal)
        throw std::runtime_error("LoadHDF5Frame frame grid size mismatch");
    ds_count.read(count.data(), H5::PredType::NATIVE_UINT8);

    file.openDataSet("vis_vx").read(vis_vx.data(), H5::PredType::NATIVE_FLOAT);
    file.openDataSet("vis_vy").read(vis_vy.data(), H5::PredType::NATIVE_FLOAT);

    file.openDataSet("vis_Jpinv").read(vis_Jpinv.data(), H5::PredType::NATIVE_FLOAT);
    file.openDataSet("vis_P").read(vis_P.data(), H5::PredType::NATIVE_FLOAT);
    file.openDataSet("vis_Q").read(vis_Q.data(), H5::PredType::NATIVE_FLOAT);

    file.openDataSet("rgb").read(rgb.data(), H5::PredType::NATIVE_UINT8);
    file.close();
    LOGR("FrameData::LoadHDF5Frame (done): {}", frameNumber);
*/
}



void FrameData::RenderFrame(VTKVisualization::VisOpt visopt)
{
    representation.ChangeVisualizationOption(visopt);

    offscreenRenderWindow->Render();
    windowToImageFilter->Modified(); // this is extra important
    std::string renderFileName = fmt::format("{:05d}.jpg", frame);

    fs::path outputDir = "output/raster";
    fs::path imageDir;
    if(visopt == VTKVisualization::VisOpt::colors) imageDir = "colors";
    else if(visopt == VTKVisualization::VisOpt::Jp_inv) imageDir = "Jpinv";
    else if(visopt == VTKVisualization::VisOpt::P) imageDir = "P";
    else if(visopt == VTKVisualization::VisOpt::Q) imageDir = "Q";
    else if(visopt == VTKVisualization::VisOpt::ridges) imageDir = "Ridges";

    fs::path targetPath = outputDir / imageDir;
    fs::create_directories(targetPath);

    fs::path fullPath = targetPath / renderFileName;

    LOGR("FrameData::RenderFrame; writing {}", renderFileName);

    writer->SetFileName(fullPath.string().c_str());
    writer->Write();
}


