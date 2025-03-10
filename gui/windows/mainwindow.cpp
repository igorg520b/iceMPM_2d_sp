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
#include "mainwindow.h"
#include "./ui_mainwindow.h"

MainWindow::~MainWindow() {delete ui;}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    params = new ParamsWrapper(&model.prms);
    representation.model = &model;
    worker = new BackgroundWorker(&model);

    // VTK
    qt_vtk_widget = new QVTKOpenGLNativeWidget();
    qt_vtk_widget->setRenderWindow(renderWindow);

    renderer->SetBackground(1.0,1.0,1.0);
    renderWindow->AddRenderer(renderer);
//    renderWindow->GetInteractor()->SetInteractorStyle(interactor);
    renderWindow->GetInteractor()->SetInteractorStyle(specialSelector2D);
    specialSelector2D->mw = this;
    renderer->AddActor(specialSelector2D->actor);

    // property browser
    pbrowser = new ObjectPropertyBrowser(this);

    // splitter
    splitter = new QSplitter(Qt::Orientation::Horizontal);
    splitter->addWidget(pbrowser);
    splitter->addWidget(qt_vtk_widget);
    splitter->setSizes(QList<int>({100, 500}));
    setCentralWidget(splitter);

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

    qsbIntentionalSlowdown = new QSpinBox();
    qsbIntentionalSlowdown->setRange(0,1000);
    qsbIntentionalSlowdown->setValue(0);
    ui->toolBar->addWidget(qsbIntentionalSlowdown);

    // statusbar
    statusLabel = new QLabel();
    labelElapsedTime = new QLabel();
    labelStepCount = new QLabel();
    labelWindSpeed = new QLabel();
    labelWindDirection = new QLabel();

    QSizePolicy sp;
    const int status_width = 90;
    sp.setHorizontalPolicy(QSizePolicy::Fixed);
    labelStepCount->setSizePolicy(sp);
    labelStepCount->setFixedWidth(status_width);
    labelElapsedTime->setSizePolicy(sp);
    labelElapsedTime->setFixedWidth(status_width);
    labelWindSpeed->setSizePolicy(sp);
    labelWindSpeed->setFixedWidth(status_width);
    labelWindDirection->setSizePolicy(sp);
    labelWindDirection->setFixedWidth(status_width);

    ui->statusbar->addWidget(statusLabel);
    ui->statusbar->addPermanentWidget(labelWindSpeed);
    ui->statusbar->addPermanentWidget(labelWindDirection);
    ui->statusbar->addPermanentWidget(labelElapsedTime);
    ui->statusbar->addPermanentWidget(labelStepCount);

// anything that includes the Model
    renderer->AddActor(representation.actor_points);
//    renderer->AddActor(representation.actor_grid);
//    renderer->AddActor(representation.actor_grid2);
    renderer->AddActor(representation.actor_uniformgrid);
    renderer->AddActor(representation.actorText);
//    renderer->AddActor(representation.scalarBar);

    // populate combobox
    QMetaEnum qme = QMetaEnum::fromType<icy::VisualRepresentation::VisOpt>();
    for(int i=0;i<qme.keyCount();i++) comboBox_visualizations->addItem(qme.key(i));

    connect(comboBox_visualizations, QOverload<int>::of(&QComboBox::currentIndexChanged),
            [&](int index){ comboboxIndexChanged_visualizations(index); });


    // slider
    slider1 = new QSlider(Qt::Horizontal);
    ui->toolBar->addWidget(slider1);
    slider1->setTracking(true);
    slider1->setMinimum(0);
    slider1->setMaximum(10000);
    connect(slider1, SIGNAL(valueChanged(int)), this, SLOT(sliderValueChanged(int)));


    // read/restore saved settings
    settingsFileName = QDir::currentPath() + "/cm.ini";
    QFileInfo fi(settingsFileName);

    if(fi.exists())
    {
        QSettings settings(settingsFileName,QSettings::IniFormat);
        QVariant var;

        vtkCamera* camera = renderer->GetActiveCamera();
        renderer->ResetCamera();
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
            memcpy(representation.ranges, ba.constData(), ba.size());
        }

        var = settings.value("take_screenshots");
        if(!var.isNull())
        {
            bool b = var.toBool();
            ui->actionTake_Screenshots->setChecked(b);
        }

        comboBox_visualizations->setCurrentIndex(settings.value("vis_option").toInt());

        var = settings.value("splitter_size_0");
        if(!var.isNull())
        {
            int sz1 = var.toInt();
            int sz2 = settings.value("splitter_size_1").toInt();
            splitter->setSizes(QList<int>({sz1, sz2}));
        }

        var = settings.value("vis_option");
        if(!var.isNull())
        {
            comboBox_visualizations->setCurrentIndex(var.toInt());
            qdsbValRange->setValue(representation.ranges[var.toInt()]);
        }
    }
    else
    {
        cameraReset_triggered();
    }


    windowToImageFilter->SetInput(renderWindow);
    windowToImageFilter->SetScale(1); // image quality
    windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
    windowToImageFilter->ReadFrontBufferOn(); // read from the back buffer
    writerPNG->SetInputConnection(windowToImageFilter->GetOutputPort());

    connect(ui->action_quit, &QAction::triggered, this, &MainWindow::quit_triggered);
    connect(ui->action_camera_reset, &QAction::triggered, this, &MainWindow::cameraReset_triggered);
    connect(ui->actionOpen, &QAction::triggered, this, &MainWindow::open_snapshot_triggered);
    connect(ui->actionStart_Pause, &QAction::triggered, this, &MainWindow::simulation_start_pause);
    connect(ui->actionLoad_Parameters, &QAction::triggered, this, &MainWindow::load_parameter_triggered);
    connect(ui->actionLoad_FLUENT_Result, &QAction::triggered, this, &MainWindow::read_fluent_data_triggered);

    connect(qdsbValRange,QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &MainWindow::limits_changed);
    connect(qsbIntentionalSlowdown,QOverload<int>::of(&QSpinBox::valueChanged), this, &MainWindow::spinbox_slowdown_value_changed);

    connect(worker, SIGNAL(workerPaused()), SLOT(background_worker_paused()));
    connect(worker, SIGNAL(stepCompleted()), SLOT(simulation_data_ready()));

    connect(params, SIGNAL(propertyChanged()), SLOT(parameters_updated()));


    pbrowser->setActiveObject(params);
    qDebug() << "MainWindow constructor done";
}


void MainWindow::closeEvent(QCloseEvent* event)
{
    quit_triggered();
    event->accept();
}


void MainWindow::quit_triggered()
{
    qDebug() << "MainWindow::quit_triggered() ";
    worker->Finalize();
    // save settings and stop simulation
    QSettings settings(settingsFileName,QSettings::IniFormat);
    qDebug() << "MainWindow: closing main window; " << settings.fileName();

    double data[10];
    renderer->GetActiveCamera()->GetPosition(&data[0]);
    renderer->GetActiveCamera()->GetFocalPoint(&data[3]);
    data[6] = renderer->GetActiveCamera()->GetParallelScale();

    qDebug() << "cam pos " << data[0] << "," << data[1] << "," << data[2];
    qDebug() << "cam focal pt " << data[3] << "," << data[4] << "," << data[5];
    qDebug() << "cam par scale " << data[6];

    QByteArray arr((char*)data, sizeof(data));
    settings.setValue("camData", arr);

    QByteArray ranges((char*)representation.ranges, sizeof(representation.ranges));
    settings.setValue("visualization_ranges", ranges);

    settings.setValue("vis_option", comboBox_visualizations->currentIndex());

    if(!qLastParameterFile.isEmpty()) settings.setValue("lastParameterFile", qLastParameterFile);

    QList<int> szs = splitter->sizes();
    settings.setValue("splitter_size_0", szs[0]);
    settings.setValue("splitter_size_1", szs[1]);

    settings.setValue("take_screenshots", ui->actionTake_Screenshots->isChecked());
    QApplication::quit();
}



void MainWindow::comboboxIndexChanged_visualizations(int index)
{
    representation.ChangeVisualizationOption(index);
    qdsbValRange->setValue(representation.ranges[index]);
    renderWindow->Render();
}

void MainWindow::limits_changed(double val_)
{
    int idx = (int)representation.VisualizingVariable;
    representation.ranges[idx] = val_;
    representation.SynchronizeValues();
    renderWindow->Render();
}

void MainWindow::cameraReset_triggered()
{
    qDebug() << "MainWindow::on_action_camera_reset_triggered()";
    vtkCamera* camera = renderer->GetActiveCamera();
    renderer->ResetCamera();
    camera->ParallelProjectionOn();
    camera->SetClippingRange(1e-1,1e3);
    camera->SetFocalPoint(0, 0., 0.);
    camera->SetPosition(0.0, 0.0, 50.0);
    camera->SetViewUp(0.0, 1.0, 0.0);
    camera->SetParallelScale(2.5);

    camera->Modified();
    renderWindow->Render();
}


void MainWindow::open_snapshot_triggered()
{
    QString defaultPath = QDir::currentPath() + "/_data/snapshots";

    // Check if the directory exists
    if (!QDir(defaultPath).exists()) {
        // Fall back to the current path if the default directory does not exist
        defaultPath = QDir::currentPath();
    }

    QString qFileName = QFileDialog::getOpenFileName(this, "Open Simulation Snapshot", defaultPath, "HDF5 Files (*.h5)");
    if(qFileName.isNull())return;
    OpenSnapshot(qFileName);
}

void MainWindow::OpenSnapshot(QString fileName)
{

    model.snapshot.ReadSnapshot(fileName.toStdString());
    representation.SynchronizeTopology();
    pbrowser->setActiveObject(params);
    updateGUI();
}


void MainWindow::load_parameter_triggered()
{
    QString qFileName = QFileDialog::getOpenFileName(this, "Load Parameters", QDir::currentPath(), "JSON Files (*.json)");
    if(qFileName.isNull())return;
    LoadParameterFile(qFileName);
}



void MainWindow::simulation_data_ready()
{
    updateGUI();
    if(ui->actionTake_Screenshots->isChecked())
        screenshot();
    model.UnlockCycleMutex();
}


void MainWindow::updateGUI()
{
    spdlog::info("updateGUI");
    labelStepCount->setText(QString::number(model.prms.SimulationStep));
    labelElapsedTime->setText(QString("%1 s").arg(model.prms.SimulationTime,0,'f',0));
//    labelWindSpeed->setText(QString("%1 m/s").arg(model.windSpeed,0,'f',2));
//    labelWindDirection->setText(QString("%1 deg").arg(model.windAngle,0,'f',0));

    // display date
    int64_t display_date = (int64_t)model.prms.SimulationTime + model.prms.SimulationStartUnixTime;
    std::time_t unix_time = display_date;
    std::tm* tm_time = std::gmtime(&unix_time);
    // Format the time
    char buffer[100];
    std::strftime(buffer, sizeof(buffer), "%d-%m-%Y %H:%M:%S UTC", tm_time);
    representation.actorText->SetInput(buffer);

    //statusLabel->setText(QString("per cycle: %1 ms").arg(model.compute_time_per_cycle,0,'f',3));

    if(model.SyncTopologyRequired)
    {
        model.SyncTopologyRequired = false;
        representation.SynchronizeTopology();
    }
    else
    {
        representation.SynchronizeValues();
    }
    renderWindow->Render();

    worker->visual_update_requested = false;
}

void MainWindow::simulation_start_pause(bool checked)
{
    if(!worker->running && checked)
    {
        qDebug() << "starting simulation via GUI";
        statusLabel->setText("starting simulation");
        worker->Resume();
    }
    else if(worker->running && !checked)
    {
        qDebug() << "pausing simulation via GUI";
        statusLabel->setText("pausing simulation");
        worker->Pause();
        ui->actionStart_Pause->setEnabled(false);
    }
}

void MainWindow::background_worker_paused()
{
    ui->actionStart_Pause->blockSignals(true);
    ui->actionStart_Pause->setEnabled(true);
    ui->actionStart_Pause->setChecked(false);
    ui->actionStart_Pause->blockSignals(false);
    statusLabel->setText("simulation stopped");
}

void MainWindow::screenshot()
{
    if(model.prms.SimulationStep % model.prms.UpdateEveryNthStep) return;
    QString outputPath = QDir::currentPath()+ "/" + screenshot_directory.c_str() + "/" +
                         QString::number(model.prms.AnimationFrameNumber()).rightJustified(5, '0') + ".png";

    QDir pngDir(QDir::currentPath()+ "/"+ screenshot_directory.c_str());
    if(!pngDir.exists()) pngDir.mkdir(QDir::currentPath()+ "/"+ screenshot_directory.c_str());

    renderWindow->DoubleBufferOff();
    renderWindow->Render();
    windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
    renderWindow->WaitForCompletion();

    windowToImageFilter->Update();
    windowToImageFilter->Modified();

    writerPNG->Modified();
    writerPNG->SetFileName(outputPath.toUtf8().constData());
    writerPNG->Write();
    renderWindow->DoubleBufferOn();
}


void MainWindow::LoadParameterFile(QString qFileName)
{
    qDebug() << "MainWindow::LoadParameterFile" << qFileName;
    std::map<std::string,std::string> additionalFiles = model.prms.ParseFile(qFileName.toStdString());

    model.snapshot.PreparePointsAndSetupGrid(additionalFiles["InputPNG"], additionalFiles["ModeledRegion"]);
    this->qLastParameterFile = qFileName;
    this->setWindowTitle(qLastParameterFile);

    if(additionalFiles.count("InputWindData")) model.snapshot.LoadWindData(additionalFiles["InputWindData"]);

    if(additionalFiles.count("InputCurrentData"))
    {
        model.fluent_interpolatror.PrepareFlowDataCache(additionalFiles["InputCurrentData"]);
    }

    if(model.prms.SaveSnapshots) {
        spdlog::info("requesting to save snapshot 0");
        model.SaveFrameRequest(model.prms.SimulationStep, model.prms.SimulationTime);
    }

    representation.SynchronizeTopology();
    pbrowser->setActiveObject(params);
    updateGUI();
}


void MainWindow::spinbox_slowdown_value_changed(int val)
{
    model.intentionalSlowdown = val;
}


void MainWindow::read_fluent_data_triggered()
{
    qDebug() << "MainWindow::read_fluent_data_triggered()";

    model.fluent_interpolatror.TestLoad(model.prms.FluentDataScale, model.prms.FluentDataOffsetX, model.prms.FluentDataOffsetY);
//    renderer->AddActor(model.fluent_interpolatror.actor);
    renderer->AddActor(model.fluent_interpolatror.actor_original);
    renderWindow->Render();



    /*
    QString defaultPath = QDir::currentPath() + "/_fluent";

    // Fall back to the current path if the default directory does not exist
    if (!QDir(defaultPath).exists()) defaultPath = QDir::currentPath();

    QString qFileName = QFileDialog::getOpenFileName(this, "Open FLUENT data", defaultPath, "FLUENT Files (*.cas.h5)");
    if(qFileName.isNull())return;
    model.fluent_interpolatror.ScanDirectory(qFileName.toStdString());

    slider1->setMinimum(0);
    slider1->setMaximum(model.fluent_interpolatror.file_count-1);

*/
}


void MainWindow::sliderValueChanged(int val)
{
    spdlog::info("sliderValueChanged {}", val);
    model.fluent_interpolatror.SetTime((double)val);
    representation.SynchronizeValues();
    renderWindow->Render();

    //model.fluent_interpolatror.LoadDataFrame(val);

 /*
    float b_val = (float)val/1000.;
    float max_val = 11.0*24*3600;
    float set_val = b_val*max_val;

    representation.wind_visualization_time = b_val*max_val;

    int64_t display_date = (int64_t)set_val + model.prms.SimulationStartUnixTime;

    std::time_t unix_time = display_date;
    std::tm* tm_time = std::gmtime(&unix_time);
    // Format the time
    char buffer[100];
    std::strftime(buffer, sizeof(buffer), "%d-%m-%Y %H:%M:%S UTC", tm_time);
    representation.actorText->SetInput(buffer);
*/
}


void MainWindow::parameters_updated()
{
    qDebug() << "MainWindow::parameters_updated(); ptsize " << model.prms.ParticleViewSize;
    representation.SynchronizeTopology();
    renderWindow->Render();
}
