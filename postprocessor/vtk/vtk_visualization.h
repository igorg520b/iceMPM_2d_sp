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

#include "framedata.h"

class VTKVisualization : public QObject
{
    Q_OBJECT

public:
    FrameData &frameData;

    VTKVisualization(FrameData& fd);

    double wind_visualization_time;

    enum VisOpt { none, count, colors, Jp_inv, velocity, P, Q, special, wind_u, wind_v, wind_norm};
    Q_ENUM(VisOpt)
    VisOpt VisualizingVariable = VisOpt::none;
    double ranges[30] = {};

    void SynchronizeValues();
    void SynchronizeTopology();
    void ChangeVisualizationOption(int option);  // invoked from GUI/main thread

    vtkNew<vtkLookupTable> hueLut, lutMPM;
    vtkNew<vtkActor> actor_grid_main;
    vtkNew<vtkActor> actor_grid_lat_lon;
    vtkNew<vtkScalarBarActor> scalarBar;
    vtkNew<vtkTextActor> actor_text;

private:
    vtkNew<vtkLookupTable> hueLut_count;
    vtkNew<vtkLookupTable> hueLut_pressure;

    // main grid
    vtkNew<vtkStructuredGrid> structuredGrid_main;
    vtkNew<vtkDataSetMapper> mapper_grid_main;
    vtkNew<vtkPoints> grid_main_points;
    vtkNew<vtkUnsignedCharArray> pts_colors;

    // lat/lon grid
    vtkNew<vtkStructuredGrid> structuredGrid_lat_lon;
    vtkNew<vtkDataSetMapper> mapper_grid_lat_lon;
    vtkNew<vtkPoints> grid_lat_lon_points;
    vtkNew<vtkFloatArray> visualized_values_grid;

    void populateLut(const float lutArray[][3], const int size, vtkNew<vtkLookupTable> &table);
    void interpolateLut(const float lutArray[][3], const int size, vtkNew<vtkLookupTable> &table);



static constexpr float lutSpecialP[10][3] = {
//    {1.,1.,1.},
    {0xb7/255.,0x24/255.,0x30/255.},
    {0xc2/255.,0x66/255.,0x6e/255.},             // -2
    {0xdf/255.,0xa7/255.,0xac/255.},             // -1

    {0xc4/255.,0xdb/255.,0xf2/255.},     // 0
    {0xc4/255.,0xdb/255.,0xf2/255.},     // 0

    {0x80/255.,0xac/255.,0xd9/255.},    //1
    {0x4c/255.,0x83/255.,0xbc/255.},    //2
    {0x27/255.,0x5d/255.,0x96/255.},    //3
    {0x0e/255.,0x41/255.,0x77/255.},    //4
    {0x03/255.,0x36/255.,0xb3/255.}
};


static constexpr float lutSpecialCount[13][3] = {
    {0x11/255.,0x18/255.,0x20/255.},    // open water
    {0.35, 0.15, 0.15},          // land
    {0.720287, 0.923781, 0.297597}, // count 1
    {0.739677, 0.775644, 0.29211},
    {0.759066, 0.627506, 0.286622},
    {0.778456, 0.479368, 0.281135},
    {0.808912, 0.363209, 0.261728},
    {0.841533, 0.316539, 0.26249},
    {0.857977, 0.321325, 0.319149},
    {0.852516, 0.319414, 0.423851},
    {0.836949, 0.282468, 0.539153},
    {0.821381, 0.245522, 0.654454},
    {0.805814, 0.208576, 0.769757}
};


static void interpolateColor(const float colorArray[][3],
                                                 int nColors, float value, float& r_out, float& g_out, float& b_out);

};
#endif
