#ifndef VTKVISUALIZATION_H
#define VTKVISUALIZATION_H

#include <QObject>

#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellType.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkNamedColors.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkLookupTable.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPoints.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkStructuredGrid.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>
#include <vtkTextActor.h>

#include <vtkRegularPolygonSource.h>
#include <vtkCylinderSource.h>

#include <vtkUnsignedCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkArrowSource.h>
#include <vtkGlyph2D.h>
#include <vtkUniformGrid.h>

#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>

#include <Eigen/Core>

#include "colormap.h"

struct FrameData;

class VTKVisualization : public QObject
{
    Q_OBJECT

public:
    FrameData *frameData;

    VTKVisualization();

    double wind_visualization_time;

    enum VisOpt { none, count, colors, Jp_inv, P, Q, wind_u, wind_v, wind_norm};
    Q_ENUM(VisOpt)
    VisOpt VisualizingVariable = VisOpt::none;
    double ranges[30] = {};

    void SynchronizeValues();
    void SynchronizeTopology();
    void ChangeVisualizationOption(int option);  // invoked from GUI/main thread

    vtkNew<vtkActor> actor_grid_main;
    vtkNew<vtkScalarBarActor> scalarBar;
    vtkNew<vtkTextActor> actor_text;
    vtkNew<vtkTextActor> actor_text_title;
    vtkNew<vtkActor> rectangleActor;

private:
    ColorMap colormap;

    vtkNew<vtkLookupTable> lut_Pressure, lut_P2;

    // main grid
    // draw background grid as vtkImageData
    vtkNew<vtkUniformGrid> uniformGrid;
    vtkNew<vtkDataSetMapper> mapper_uniformgrid;

    // frame around the grid
    vtkNew<vtkPolyLine> polyline;
    vtkNew<vtkPolyData> rectanglePolyData;
    vtkNew<vtkPoints> rectanglePoints;
    vtkNew<vtkCellArray> rectangleLines;
    vtkNew<vtkPolyDataMapper> rectangleMapper;

    vtkNew<vtkWindowToImageFilter> windowToImageFilter;
    vtkNew<vtkPNGWriter> writerPNG;
    vtkNew<vtkRenderer> offscreenRenderer;

};
#endif
