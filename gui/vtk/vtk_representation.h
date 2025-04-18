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
#include <vtkPlaneSource.h>
#include <vtkTexture.h>

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
    vtkNew<vtkActor> raster_actor;
    vtkNew<vtkActor> testing_actor;

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

    // background image
    std::vector<uint8_t> renderedImage;
    vtkNew<vtkImageData> raster_imageData;
    vtkNew<vtkUnsignedCharArray> raster_scalars;
    vtkNew<vtkPlaneSource> raster_plane;
    vtkNew<vtkTexture> raster_texture;
    vtkNew<vtkPolyDataMapper> raster_mapper;


    // one point per cell (in the center) for testing
    vtkNew<vtkPoints> testing_points;
    vtkNew<vtkPolyData> testing_points_polydata;
    vtkNew<vtkPolyDataMapper> testing_points_mapper;
    vtkNew<vtkVertexGlyphFilter> testing_points_filter;


    static constexpr uint8_t waterColor[3] = {0x15, 0x1f, 0x2f};
};
#endif
