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
    : QMainWindow(parent), frameData(ggd)
    , ui(new Ui::PPMainWindow)
{
    ui->setupUi(this);

    // VTK
    qt_vtk_widget = new QVTKOpenGLNativeWidget();
    qt_vtk_widget->setRenderWindow(renderWindow);
    setCentralWidget(qt_vtk_widget);


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


    slider2 = new QSlider(Qt::Horizontal);
    ui->toolBar->addWidget(slider2);
    slider2->setTracking(true);
    slider2->setMinimum(1);
    slider2->setMaximum(10000);
    connect(slider2, SIGNAL(valueChanged(int)), this, SLOT(sliderValueChanged(int)));

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



    connect(ui->action_camera_reset, &QAction::triggered, this, &PPMainWindow::cameraReset_triggered);

//    connect(ui->actionOpen_Frame, &QAction::triggered, this, &PPMainWindow::open_frame_triggered);
//    connect(ui->actionRender_Frame, &QAction::triggered, this, &PPMainWindow::render_frame_triggered);
//    connect(ui->actionGenerate_Script, &QAction::triggered, this, &PPMainWindow::generate_script_triggered);
//    connect(ui->actionShow_Wind, &QAction::triggered, this, &PPMainWindow::toggle_wind_visualization);
//    connect(ui->actionRender_All, &QAction::triggered, this, &PPMainWindow::render_all_triggered);
//    connect(ui->actionRender_All_2, &QAction::triggered, this, &PPMainWindow::render_all2_triggered);
//    connect(ui->actionLoad_Selected_Frame, &QAction::triggered, this, &PPMainWindow::load_selected_frame_triggered);

    connect(qdsbValRange,QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &PPMainWindow::limits_changed);

    // off-screen rendering
    // Copy the renderers from the current render window

    ui->actionShow_Wind->setChecked(false);
    // Set the target resolution for the off-screen render window

    qDebug() << "PPMainWindow constructor done";
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




// ==================================================

void PPMainWindow::LoadParametersFile(QString fileName)
{
    qDebug() << "PPMainWindow::LoadParametersFile: " << fileName;
    ggd.ReadParameterFile(fileName.toStdString());
}

void PPMainWindow::LoadFramesDirectory(QString framesDirectory)
{
    qDebug() << "PPMainWindow::LoadFramesDirectory: " << framesDirectory;
    ggd.ScanDirectory(framesDirectory.toStdString());
    slider2->setMaximum(ggd.countFrames);
    frameData.LoadHDF5Frame(1);
    frameData.representation.SynchronizeValues();
    renderWindow->Render();
}

void PPMainWindow::sliderValueChanged(int val)
{
    qDebug() << "slider set to " << val;
    frameData.LoadHDF5Frame(val);
    frameData.representation.SynchronizeValues();
    renderWindow->Render();

    renderWindow->Render();
}


void PPMainWindow::comboboxIndexChanged_visualizations(int index)
{
    frameData.representation.ChangeVisualizationOption(index);
    qdsbValRange->setValue(frameData.representation.ranges[index]);
    renderWindow->Render();
}


/*
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

    qsbFrameFrom->setRange(1, ggd.countFrames-1);
    qsbFrameFrom->setValue(1);
    qsbFrameTo->setRange(1, ggd.countFrames-1);
    qsbFrameTo->setValue(ggd.countFrames-1);
    qsbLoadFrame->setRange(1, ggd.countFrames-1);

    updateGUI();
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
    localFD.RenderFrame(frameData.representation.VisualizingVariable,
                        getOutputDirectory(frameData.representation.VisualizingVariable));
}


void PPMainWindow::render_all_triggered()
{
    qDebug() << "PPMainWindow::render_all_triggered()";
    const int frameFrom = qsbFrameFrom->value();
    const int frameTo = qsbFrameTo->value();

    std::mutex mutex;
    int threads = qsbThreads->value();

//    omp_set_num_threads(threads);

    double data[10];
    frameData.renderer->GetActiveCamera()->GetPosition(&data[0]);
    frameData.renderer->GetActiveCamera()->GetFocalPoint(&data[3]);
    data[6] = frameData.renderer->GetActiveCamera()->GetParallelScale();

    FrameData localFD;
    localFD.mutex = &mutex;

//
//#pragma omp parallel for schedule(dynamic, 5) firstprivate(localFD) num_threads(2)
    for(int frame=frameFrom; frame<=frameTo; frame++)
    {
        localFD.SetUpOffscreenRender(frameData, data);
        localFD.LoadHDF5Frame(frame);
        localFD.RenderFrame(frameData.representation.VisualizingVariable,
                            getOutputDirectory(frameData.representation.VisualizingVariable));
    }
}


void PPMainWindow::render_all2_triggered()
{
    qDebug() << "PPMainWindow::render_all_triggered()";
    const int frameFrom = qsbFrameFrom->value();
    const int frameTo = qsbFrameTo->value();

    std::mutex mutex;
    int threads = qsbThreads->value();

    //    omp_set_num_threads(threads);

    double data[10];
    frameData.renderer->GetActiveCamera()->GetPosition(&data[0]);
    frameData.renderer->GetActiveCamera()->GetFocalPoint(&data[3]);
    data[6] = frameData.renderer->GetActiveCamera()->GetParallelScale();

    FrameData localFD;
    localFD.mutex = &mutex;

//
#pragma omp parallel for schedule(dynamic, 5) firstprivate(localFD) num_threads(10)
    for(int frame=frameFrom; frame<=frameTo; frame++)
    {
        localFD.SetUpOffscreenRender(frameData, data);
        localFD.LoadHDF5Frame(frame);
        localFD.RenderFrame(VTKVisualization::VisOpt::Jp_inv,
                            getOutputDirectory(VTKVisualization::VisOpt::Jp_inv));
        localFD.RenderFrame(VTKVisualization::VisOpt::colors,
                            getOutputDirectory(VTKVisualization::VisOpt::colors));
        localFD.RenderFrame(VTKVisualization::VisOpt::P,
                            getOutputDirectory(VTKVisualization::VisOpt::P));
        localFD.RenderFrame(VTKVisualization::VisOpt::Q,
                            getOutputDirectory(VTKVisualization::VisOpt::Q));
    }
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

*/
