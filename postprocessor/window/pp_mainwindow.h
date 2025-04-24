#ifndef PPMAINWINDOW_H
#define PPMAINWINDOW_H

#include <QMainWindow>

#include <QSizePolicy>
#include <QPushButton>
#include <QSplitter>
#include <QLabel>
#include <QVBoxLayout>
#include <QTreeWidget>
#include <QProgressBar>
#include <QMenu>
#include <QList>
#include <QDebug>
#include <QComboBox>
#include <QMetaEnum>
#include <QDir>
#include <QString>
#include <QCheckBox>
#include <QFile>
#include <QTextStream>
#include <QIODevice>
#include <QSettings>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QFileInfo>

#include <QVTKOpenGLNativeWidget.h>

#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkProperty.h>
#include <vtkNew.h>
#include <vtkInteractorStyleImage.h>

#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>

#include "vtk_visualization.h"
#include "framedata.h"


QT_BEGIN_NAMESPACE
namespace Ui { class PPMainWindow; }
QT_END_NAMESPACE

class PPMainWindow : public QMainWindow
{
    Q_OBJECT
private:
    Ui::PPMainWindow *ui;

public:
    PPMainWindow(QWidget *parent = nullptr);
    ~PPMainWindow();

    void closeEvent( QCloseEvent* event ) override;

private Q_SLOTS:
    void cameraReset_triggered();
    void comboboxIndexChanged_visualizations(int index);
    void limits_changed(double val);

    void open_frame_triggered();
    void render_frame_triggered();
    void render_all_triggered();
    void render_all2_triggered();
    void generate_script_triggered();
    void load_selected_frame_triggered();

    void toggle_wind_visualization(bool checked);

private:
    FrameData frameData;
    GeneralGridData ggd;


    void updateGUI();
    void OpenSnapshot(QString fileName);

    QString settingsFileName;       // includes current dir
    QComboBox *comboBox_visualizations;
    QDoubleSpinBox *qdsbValRange;   // high and low limits for value scale
    QSpinBox *qsbFrameFrom, *qsbFrameTo, *qsbThreads, *qsbLoadFrame;
    QLabel *lbl;

    // VTK
    vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
    QVTKOpenGLNativeWidget *qt_vtk_widget;

    // other
    vtkNew<vtkInteractorStyleImage> interactor;

    constexpr static std::string_view outputDirectoryP = "render/P";
    constexpr static std::string_view outputDirectoryQ = "render/Q";
    constexpr static std::string_view outputDirectoryColors = "render/colors";
    constexpr static std::string_view outputDirectoryJpinv = "render/Jp_inv";

    static std::string_view getOutputDirectory(VTKVisualization::VisOpt visopt);
};
#endif
