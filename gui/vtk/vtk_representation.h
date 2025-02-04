#ifndef MESHREPRESENTATION_H
#define MESHREPRESENTATION_H

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
#include <vtkImageData.h>
#include <vtkUniformGrid.h>

#include "colormap.h"

namespace icy { class VisualRepresentation; class Model;}

class icy::VisualRepresentation : public QObject
{
    Q_OBJECT

public:
    VisualRepresentation();

    icy::Model *model;
    double wind_visualization_time;

    enum VisOpt { none, status, Jp_inv, P, Q, color, v_u, v_v, v_norm, strength};
    Q_ENUM(VisOpt)
    VisOpt VisualizingVariable = VisOpt::none;
    double ranges[30] = {};

    void SynchronizeValues();
    void SynchronizeTopology();
    void ChangeVisualizationOption(int option);  // invoked from GUI/main thread

    vtkNew<vtkActor> actor_points;
    vtkNew<vtkActor> actor_grid;
    vtkNew<vtkActor> actor_uniformgrid;
    vtkNew<vtkScalarBarActor> scalarBar;
    vtkNew<vtkTextActor> actorText;

private:
    ColorMap colormap;

    // points
    vtkNew<vtkPoints> points;
    vtkNew<vtkPolyData> points_polydata;
    vtkNew<vtkPolyDataMapper> points_mapper;

    vtkNew<vtkVertexGlyphFilter> points_filter;
    vtkNew<vtkUnsignedCharArray> pts_colors;    // for color visualization per point

    // draw background grid as vtkImageData
    vtkNew<vtkUniformGrid> uniformGrid;
    vtkNew<vtkDataSetMapper> mapper_uniformgrid;

    // background grid
    vtkNew<vtkStructuredGrid> structuredGrid;
    vtkNew<vtkDataSetMapper> grid_mapper;
    vtkNew<vtkPoints> grid_points;
};
#endif
