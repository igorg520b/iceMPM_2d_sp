#include <filesystem>
#include <regex>
#include "framedata.h"




void GeneralGridData::ScanDirectory(std::string frameFileName)
{
    std::filesystem::path framePath(frameFileName);

    // Navigate to the target file by appending the relative path
    std::filesystem::path parentPath = framePath.parent_path();
    frameDirectory = parentPath.string();

    std::vector<bool> availableFrames;

    // Read the contents of the directory
    for (const auto& entry : std::filesystem::directory_iterator(frameDirectory)) {
        if (entry.is_regular_file()) {
            std::string filename = entry.path().filename().string();

            // Check if the filename matches the pattern "frame_XXXXX.h5"
            std::regex pattern(R"(frame_(\d{5})\.h5)");
            std::smatch match;
            if (std::regex_match(filename, match, pattern)) {
                // Extract the five-digit integer from the filename
                int frameIndex = std::stoi(match[1].str());

                // Ensure the availableFrames vector is large enough
                if (frameIndex >= availableFrames.size()) {
                    availableFrames.resize(frameIndex + 1, false);
                }

                // Set the corresponding element to true
                availableFrames[frameIndex] = true;
            }
        }
    }

    countFrames = availableFrames.size();
    LOGR("availableFrames size {}",availableFrames.size());
    LOGR("frameDirectory {}", frameDirectory);

    // Load grid data

    // Navigate to the target file by appending the relative path
    std::filesystem::path gridAndWindPath = parentPath / "../snapshots/_grid_and_wind.h5";

    // Normalize the path (resolving ".." and ".")
    gridAndWindPath = std::filesystem::canonical(gridAndWindPath);

    // Convert back to std::string if needed
    std::string resultPath = gridAndWindPath.string();

    // open the initial data file
    LOGR("LoadHDF5Frame resultPath = {}", resultPath);

    H5::H5File file(resultPath, H5F_ACC_RDONLY);

    H5::DataSet gridDataset = file.openDataSet("grid");
    prms.ReadParametersFromHDF5Attributes(gridDataset);
    prms.Printout();

    // Get the dataspace of the dataset
    H5::DataSpace dataspace = gridDataset.getSpace();

    // Get the dimensions of the dataset
    hsize_t dims[2];
    dataspace.getSimpleExtentDims(dims, nullptr);
    size_t gridXTotal = dims[0];
    size_t gridY = dims[1];

    if(gridXTotal != prms.GridXTotal || gridY != prms.GridYTotal)
        throw std::runtime_error("LoadHDF5Frame grid size mismatch");

    // Resize the vector to hold the data
    grid_status_buffer.resize(gridXTotal * gridY);

    // Read the dataset into the vector
    gridDataset.read(grid_status_buffer.data(), H5::PredType::NATIVE_UINT8);

    // read original image color
    original_image_colors_rgb.resize(prms.InitializationImageSizeX*prms.InitializationImageSizeY*3);
    H5::DataSet rgbDataset = file.openDataSet("colorImage");
    rgbDataset.read(original_image_colors_rgb.data(), H5::PredType::NATIVE_UINT8);

    LOGV("GeneralGridData::ScanDirectory done");
}




//==========================================================================

FrameData::FrameData(const FrameData &other)
{
    this->mutex = other.mutex;
    representation.frameData = this;

    offscreenRenderWindow->Initialize();

    renderer->SetBackground(1.0,1.0,1.0);
    offscreenRenderWindow->SetSize(1920, 1080);
    offscreenRenderWindow->DoubleBufferOff();
    offscreenRenderWindow->SetOffScreenRendering(true);

    renderer->AddActor(representation.actor_grid_main);
    renderer->AddActor(representation.actor_text);
    renderer->AddActor(representation.scalarBar);
    renderer->AddActor(representation.rectangleActor);
    renderer->AddActor(representation.actor_text_title);

    windowToImageFilter->SetInput(offscreenRenderWindow);
    windowToImageFilter->SetScale(1); // image quality
    windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
    windowToImageFilter->ReadFrontBufferOn(); // read from the back buffer
    writerPNG->SetInputConnection(windowToImageFilter->GetOutputPort());
}



FrameData::FrameData()
{
    representation.frameData = this;

    offscreenRenderWindow->Initialize();

    renderer->SetBackground(1.0,1.0,1.0);
    offscreenRenderWindow->SetSize(1920, 1080);
    offscreenRenderWindow->DoubleBufferOff();
    offscreenRenderWindow->SetOffScreenRendering(true);

    renderer->AddActor(representation.actor_grid_main);
    renderer->AddActor(representation.actor_text);
    renderer->AddActor(representation.scalarBar);
    renderer->AddActor(representation.rectangleActor);
    renderer->AddActor(representation.actor_text_title);

    windowToImageFilter->SetInput(offscreenRenderWindow);
    windowToImageFilter->SetScale(1); // image quality
    windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
    windowToImageFilter->ReadFrontBufferOn(); // read from the back buffer
    writerPNG->SetInputConnection(windowToImageFilter->GetOutputPort());
}


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


void FrameData::LoadHDF5Frame(int frameNumber)
{
    LOGR("LoadHDF5Frame: {}; ggd->frameDirectory {}", frameNumber, ggd->frameDirectory);
    std::string fileName = fmt::format("{}/frame_{:05d}.h5", ggd->frameDirectory, frameNumber);
    mutex->lock();
    LoadHDF5Frame(fileName);
    mutex->unlock();
    this->frame = frameNumber;
}


void FrameData::LoadHDF5Frame(std::string fileName)
{
    LOGR("invoked LoadHDF5Frame(std::string fileName) {}", fileName);
    size_t gridSize = ggd->prms.GridXTotal*ggd->prms.GridYTotal;
    LOGR("LoadHDF5Frame gridSize {}x{} = {}",ggd->prms.GridXTotal, ggd->prms.GridYTotal, gridSize );
    count.resize(gridSize);
    vis_vx.resize(gridSize);
    vis_vy.resize(gridSize);
    vis_Jpinv.resize(gridSize);
    vis_P.resize(gridSize);
    vis_Q.resize(gridSize);
    rgb.resize(gridSize*3);

    H5::H5File file(fileName, H5F_ACC_RDONLY);
    H5::DataSet ds_count = file.openDataSet("count");

    ds_count.openAttribute("SimulationStep").read(H5::PredType::NATIVE_INT, &SimulationStep);
    ds_count.openAttribute("SimulationTime").read(H5::PredType::NATIVE_DOUBLE, &SimulationTime);

    hsize_t dims[2];
    ds_count.getSpace().getSimpleExtentDims(dims, nullptr);
    if(dims[0] != ggd->prms.GridXTotal || dims[1] != ggd->prms.GridYTotal)
        throw std::runtime_error("LoadHDF5Frame frame grid size mismatch");
    ds_count.read(count.data(), H5::PredType::NATIVE_UINT8);

    file.openDataSet("vis_vx").read(vis_vx.data(), H5::PredType::NATIVE_FLOAT);
    file.openDataSet("vis_vy").read(vis_vy.data(), H5::PredType::NATIVE_FLOAT);

    file.openDataSet("vis_Jpinv").read(vis_Jpinv.data(), H5::PredType::NATIVE_FLOAT);
    file.openDataSet("vis_P").read(vis_P.data(), H5::PredType::NATIVE_FLOAT);
    file.openDataSet("vis_Q").read(vis_Q.data(), H5::PredType::NATIVE_FLOAT);

    file.openDataSet("rgb").read(rgb.data(), H5::PredType::NATIVE_UINT8);
    dataLoaded = true;
    LOGV("LoadHDF5Frame done");
}


