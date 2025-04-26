#include "framedata.h"
#include <fmt/format.h>
#include <fmt/std.h>



FrameData::FrameData(GeneralGridData &ggd_) : ggd(ggd_), representation(*this)
{
//    offscreenRenderWindow->Initialize();

    //renderer->SetBackground(1.0,1.0,1.0);
//    offscreenRenderWindow->SetSize(1920, 1080);
//    offscreenRenderWindow->DoubleBufferOff();
//    offscreenRenderWindow->SetOffScreenRendering(true);

//    renderer->AddActor(representation.raster_actor);
//    renderer->AddActor(representation.actor_text);
//    renderer->AddActor(representation.scalarBar);
//    renderer->AddActor(representation.actor_text_title);

//    windowToImageFilter->SetInput(offscreenRenderWindow);
//    windowToImageFilter->SetScale(1); // image quality
//    windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
//    windowToImageFilter->ReadFrontBufferOn(); // read from the back buffer
//    writerPNG->SetInputConnection(windowToImageFilter->GetOutputPort());
}


void FrameData::LoadHDF5Frame(int frameNumber)
{
    this->frame = frameNumber;
    //f00001.h5
    std::string fileName = fmt::format(fmt::runtime("{}/f{:05d}.h5"), ggd.frameDirectory, frame);

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
}




/*

void FrameData::SetUpOffscreenRender(FrameData &guiFD, double data[10])
{
    this->ggd = guiFD.ggd;

    // copy ranges array
    const VTKVisualization& other = guiFD.representation;
    std::copy(std::begin(other.ranges), std::end(other.ranges), std::begin(this->representation.ranges));

    representation.VisualizingVariable = guiFD.representation.VisualizingVariable;

    // set up camera

    vtkCamera* camera = renderer->GetActiveCamera();
    renderer->ResetCamera();
    camera->ParallelProjectionOn();

    camera->SetClippingRange(1e-1,1e4);
    camera->SetViewUp(0.0, 1.0, 0.0);
    camera->SetPosition(data[0],data[1],data[2]);
    camera->SetFocalPoint(data[3],data[4],data[5]);
    camera->SetParallelScale(data[6]);
    camera->Modified();

    // renderer -> offscreen
    offscreenRenderWindow->AddRenderer(renderer);
}

void FrameData::RenderFrame(VTKVisualization::VisOpt visopt, std::string_view outputDirectory)
{
    LOGV("FrameData::RenderFrame");
    representation.ChangeVisualizationOption(visopt);
    offscreenRenderWindow->Render();
    windowToImageFilter->Modified(); // this is extra important
    std::string renderFileName = fmt::format("{}/{:05d}.png", outputDirectory, frame);
    LOGR("writing {}", renderFileName);
    writerPNG->SetFileName(renderFileName.c_str());
    mutex->lock();
    writerPNG->Write();
    mutex->unlock();
}


void FrameData::LoadHDF5Frame(std::string fileName)
{

}

*/
