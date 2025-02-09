#include <QFileDialog>
#include <QList>
#include <QPointF>
#include <QCloseEvent>
#include <QStringList>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <filesystem>

#include <spdlog/spdlog.h>
#include <omp.h>

#include "pp_mainwindow.h"
#include "./ui_pp_mainwindow.h"



PPMainWindow::~PPMainWindow() {delete ui;}

PPMainWindow::PPMainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::PPMainWindow)
{
    ui->setupUi(this);

    // VTK
    qt_vtk_widget = new QVTKOpenGLNativeWidget();
    qt_vtk_widget->setRenderWindow(renderWindow);
    setCentralWidget(qt_vtk_widget);


    frameData.ggd = &ggd;
    renderWindow->AddRenderer(frameData.renderer);
    renderWindow->GetInteractor()->SetInteractorStyle(interactor);

    // toolbar - combobox
    comboBox_visualizations = new QComboBox();
    ui->toolBar->addWidget(comboBox_visualizations);

    // double spin box
    qdsbValRange = new QDoubleSpinBox();
    qdsbValRange->setRange(-10, 10);
    qdsbValRange->setValue(-2);
    qdsbValRange->setDecimals(2);
    qdsbValRange->setSingleStep(0.25);
    ui->toolBar->addWidget(qdsbValRange);

    qsbFrameFrom = new QSpinBox();
    qsbFrameTo = new QSpinBox();
    ui->toolBar->addWidget(qsbFrameFrom);
    ui->toolBar->addWidget(qsbFrameTo);

    // populate combobox
    QMetaEnum qme = QMetaEnum::fromType<VTKVisualization::VisOpt>();
    for(int i=0;i<qme.keyCount();i++) comboBox_visualizations->addItem(qme.key(i));

    connect(comboBox_visualizations, QOverload<int>::of(&QComboBox::currentIndexChanged),
            [&](int index){ comboboxIndexChanged_visualizations(index); });



    // read/restore saved settings
    settingsFileName = QDir::currentPath() + "/ppcm.ini";
    QFileInfo fi(settingsFileName);

    if(fi.exists())
    {
        QSettings settings(settingsFileName,QSettings::IniFormat);
        QVariant var;

        vtkCamera* camera = frameData.renderer->GetActiveCamera();
        frameData.renderer->ResetCamera();
        camera->ParallelProjectionOn();

        var = settings.value("camData");
        if(!var.isNull())
        {
            double *vec = (double*)var.toByteArray().constData();
            camera->SetClippingRange(1e-1,1e4);
            camera->SetViewUp(0.0, 1.0, 0.0);
            camera->SetPosition(vec[0],vec[1],vec[2]);
            camera->SetFocalPoint(vec[3],vec[4],vec[5]);
            camera->SetParallelScale(vec[6]);
            camera->Modified();
        }

        var = settings.value("visualization_ranges");
        if(!var.isNull())
        {
            QByteArray ba = var.toByteArray();
            memcpy(frameData.representation.ranges, ba.constData(), ba.size());
        }

        var = settings.value("vis_option");
        if(!var.isNull())
        {
            comboBox_visualizations->setCurrentIndex(var.toInt());
            qdsbValRange->setValue(frameData.representation.ranges[var.toInt()]);
        }
    }
    else
    {
        cameraReset_triggered();
    }

    connect(ui->actionOpen_Frame, &QAction::triggered, this, &PPMainWindow::open_frame_triggered);
    connect(ui->action_camera_reset, &QAction::triggered, this, &PPMainWindow::cameraReset_triggered);
    connect(ui->actionRender_Frame, &QAction::triggered, this, &PPMainWindow::render_frame_triggered);
    connect(ui->actionRender_All, &QAction::triggered, this, &PPMainWindow::render_all_triggered);
    connect(ui->actionGenerate_Script, &QAction::triggered, this, &PPMainWindow::generate_script_triggered);
    connect(ui->actionShow_Wind, &QAction::triggered, this, &PPMainWindow::toggle_wind_visualization);

    connect(qdsbValRange,QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &PPMainWindow::limits_changed);

    // off-screen rendering
    // Copy the renderers from the current render window

    ui->actionShow_Wind->setChecked(false);
    // Set the target resolution for the off-screen render window

    std::filesystem::path directory_path1(outputDirectoryP);
    if (!std::filesystem::exists(directory_path1)) std::filesystem::create_directories(directory_path1);
    std::filesystem::path directory_path2(outputDirectoryQ);
    if (!std::filesystem::exists(directory_path2)) std::filesystem::create_directories(directory_path2);
    std::filesystem::path directory_path3(outputDirectoryColors);
    if (!std::filesystem::exists(directory_path3)) std::filesystem::create_directories(directory_path3);
    std::filesystem::path directory_path4(outputDirectoryJpinv);
    if (!std::filesystem::exists(directory_path4)) std::filesystem::create_directories(directory_path4);


    qDebug() << "PPMainWindow constructor done";
    updateGUI();
}

void PPMainWindow::updateGUI()
{

    renderWindow->Render();
}

void PPMainWindow::open_frame_triggered()
{
    QString defaultPath = QDir::currentPath() + "/_data/frames";

    // Check if the directory exists
    if (!QDir(defaultPath).exists()) {
        // Fall back to the current path if the default directory does not exist
        defaultPath = QDir::currentPath();
    }

    // Open the file dialog
    QString qFileName = QFileDialog::getOpenFileName(
        this,
        "Open Simulation Snapshot",
        defaultPath,
        "HDF5 Files (*.h5)"
        );

    if(qFileName.isNull())return;
    ggd.ScanDirectory(qFileName.toStdString());


    frameData.LoadHDF5Frame(qFileName.toStdString());
    frameData.representation.SynchronizeTopology();

    qsbFrameFrom->setRange(0, ggd.countFrames-1);
    qsbFrameFrom->setValue(0);
    qsbFrameTo->setRange(0, ggd.countFrames-1);
    qsbFrameTo->setValue(ggd.countFrames-1);

    updateGUI();
}


void PPMainWindow::comboboxIndexChanged_visualizations(int index)
{
    frameData.representation.ChangeVisualizationOption(index);
    qdsbValRange->setValue(frameData.representation.ranges[index]);
    renderWindow->Render();
}

void PPMainWindow::closeEvent(QCloseEvent* event)
{
    qDebug() << "close event";

    QSettings settings(settingsFileName,QSettings::IniFormat);
    qDebug() << "PPMainWindow: closing main window; " << settings.fileName();

    double data[10];
    frameData.renderer->GetActiveCamera()->GetPosition(&data[0]);
    frameData.renderer->GetActiveCamera()->GetFocalPoint(&data[3]);
    data[6] = frameData.renderer->GetActiveCamera()->GetParallelScale();

    qDebug() << "cam pos " << data[0] << "," << data[1] << "," << data[2];
    qDebug() << "cam focal pt " << data[3] << "," << data[4] << "," << data[5];
    qDebug() << "cam par scale " << data[6];

    QByteArray arr((char*)data, sizeof(data));
    settings.setValue("camData", arr);

    QByteArray ranges((char*)frameData.representation.ranges, sizeof(frameData.representation.ranges));
    settings.setValue("visualization_ranges", ranges);

    settings.setValue("vis_option", comboBox_visualizations->currentIndex());

    event->accept();
}


void PPMainWindow::limits_changed(double val_)
{
    qDebug() << "limits_changed";
    int idx = (int)frameData.representation.VisualizingVariable;
    frameData.representation.ranges[idx] = val_;
    frameData.representation.SynchronizeValues();
    renderWindow->Render();
}



void PPMainWindow::cameraReset_triggered()
{
    qDebug() << "MainWindow::on_action_camera_reset_triggered()";
    vtkCamera* camera = frameData.renderer->GetActiveCamera();
    frameData.renderer->ResetCamera();
    camera->ParallelProjectionOn();
    camera->SetClippingRange(1e-1,1e3);

    const double dx = ggd.prms.cellsize*ggd.prms.InitializationImageSizeX/2;
    const double dy = ggd.prms.cellsize*ggd.prms.InitializationImageSizeY/2;

    qDebug() << "dx " << dx << "\ndy " << dy;

    camera->SetPosition(dx, dy, 50.);
    camera->SetFocalPoint(dx, dy, 0.);
    camera->SetViewUp(0.0, 1.0, 0.0);
    camera->SetParallelScale(std::min(dx,dy)*1.1);

    camera->Modified();
    renderWindow->Render();
}


void PPMainWindow::render_frame_triggered()
{
    qDebug() << "PPMainWindow::render_frame_triggered()";
    FrameData localFD;

    double data[10];
    frameData.renderer->GetActiveCamera()->GetPosition(&data[0]);
    frameData.renderer->GetActiveCamera()->GetFocalPoint(&data[3]);
    data[6] = frameData.renderer->GetActiveCamera()->GetParallelScale();


    localFD.SetUpOffscreenRender(frameData, data);
    localFD.LoadHDF5Frame(ggd.countFrames/2);
    localFD.RenderFrame(VTKVisualization::VisOpt::Jp_inv, outputDirectoryJpinv);
}



void PPMainWindow::render_all_triggered()
{
    qDebug() << "PPMainWindow::render_all_triggered()";
    const int frameFrom = qsbFrameFrom->value();
    const int frameTo = qsbFrameTo->value();

    std::mutex mutex;

//    omp_set_num_threads(6);

    double data[10];
    frameData.renderer->GetActiveCamera()->GetPosition(&data[0]);
    frameData.renderer->GetActiveCamera()->GetFocalPoint(&data[3]);
    data[6] = frameData.renderer->GetActiveCamera()->GetParallelScale();


#pragma omp parallel for schedule(dynamic, 1) num_threads(1)
    for(int frame=frameFrom; frame<=frameTo; frame++)
    {
        FrameData localFD;
        localFD.mutex = &mutex;
        localFD.SetUpOffscreenRender(frameData, data);
        localFD.LoadHDF5Frame(frame);
        localFD.RenderFrame(VTKVisualization::VisOpt::Jp_inv, outputDirectoryJpinv);
    }


    /*


    renderWindow->RemoveRenderer(renderer);
    offscreenRenderWindow->AddRenderer(renderer);






        if(!frameData.availableFrames[frame]) continue;
        qDebug() << "rendering frame " << frame;
        frameData.LoadHDF5Frame(fileName, false);

        std::string renderFileName;

        representation.ChangeVisualizationOption(VTKVisualization::VisOpt::Jp_inv);
        offscreenRenderWindow->Render();
        windowToImageFilter->Modified(); // this is extra important
        renderFileName = fmt::format("{}/{:05d}.png", outputDirectoryJpinv, frame);
        writerPNG->SetFileName(renderFileName.c_str());
        writerPNG->Write();


        representation.ChangeVisualizationOption(VTKVisualization::VisOpt::P);
        offscreenRenderWindow->Render();
        windowToImageFilter->Modified(); // this is extra important
        std::string renderFileName = fmt::format("{}/{:05d}.png", outputDirectoryP, frame);
        writerPNG->SetFileName(renderFileName.c_str());
        writerPNG->Write();

        representation.ChangeVisualizationOption(VTKVisualization::VisOpt::Q);
        offscreenRenderWindow->Render();
        windowToImageFilter->Modified(); // this is extra important
        renderFileName = fmt::format("{}/{:05d}.png", outputDirectoryQ, frame);
        writerPNG->SetFileName(renderFileName.c_str());
        writerPNG->Write();


        representation.ChangeVisualizationOption(VTKVisualization::VisOpt::colors);
        offscreenRenderWindow->Render();
        windowToImageFilter->Modified(); // this is extra important
        renderFileName = fmt::format("{}/{:05d}.png", outputDirectoryColors, frame);
        writerPNG->SetFileName(renderFileName.c_str());
        writerPNG->Write();

    }

    offscreenRenderWindow->RemoveRenderer(renderer);
    renderWindow->AddRenderer(renderer);
*/
}


void PPMainWindow::generate_script_triggered()
{
    std::string filename = "render/genvideo.sh";
    std::ofstream scriptFile(filename);
    int frames = ggd.countFrames;

    std::string cmd = fmt::format("ffmpeg -y -r 30 -f image2 -start_number 0 -i \"P/%05d.png\" -vframes {} -vcodec libx264 -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -crf 21  -pix_fmt yuv420p \"P.mp4\"", frames-1);
    scriptFile << cmd << '\n';
    cmd = fmt::format("ffmpeg -y -r 30 -f image2 -start_number 0 -i \"Q/%05d.png\" -vframes {} -vcodec libx264 -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -crf 21  -pix_fmt yuv420p \"Q.mp4\"", frames-1);
    scriptFile << cmd << '\n';
    cmd = fmt::format("ffmpeg -y -r 30 -f image2 -start_number 0 -i \"Jp_inv/%05d.png\" -vframes {} -vcodec libx264 -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -crf 21  -pix_fmt yuv420p \"Jp_inv.mp4\"", frames-1);
    scriptFile << cmd << '\n';
    cmd = fmt::format("ffmpeg -y -r 30 -f image2 -start_number 0 -i \"colors/%05d.png\" -vframes {} -vcodec libx264 -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -crf 21  -pix_fmt yuv420p \"colors.mp4\"", frames-1);
    scriptFile << cmd << '\n';

    std::string concat = R"(ffmpeg -i P.mp4 -i Q.mp4 -i Jp_inv.mp4 -i colors.mp4 -filter_complex "[0:v:0][1:v:0][2:v:0][3:v:0]concat=n=4:v=1[outv]" -map "[outv]" output.mp4)";
    scriptFile << concat;
    scriptFile.close();

    // Set the script as executable
    if (std::system(("chmod +x " + filename).c_str()) != 0) {
        std::cerr << "Error: Failed to set executable permission\n";
    }
}


void PPMainWindow::toggle_wind_visualization(bool checked)
{
    updateGUI();
}
