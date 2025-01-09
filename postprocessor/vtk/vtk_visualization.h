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


#include <Eigen/Core>

#include "framedata.h"

class VTKVisualization : public QObject
{
    Q_OBJECT

public:
    FrameData &frameData;

    VTKVisualization(FrameData& fd);

    double wind_visualization_time;

    enum VisOpt { none, count, colors, Jp_inv, P, Q, wind_u, wind_v, wind_norm};
    Q_ENUM(VisOpt)
    VisOpt VisualizingVariable = VisOpt::none;
    double ranges[30] = {};

    void SynchronizeValues();
    void SynchronizeTopology();
    void ChangeVisualizationOption(int option);  // invoked from GUI/main thread

    vtkNew<vtkActor> actor_grid_lat_lon;

    vtkNew<vtkActor> actor_grid_main;
    vtkNew<vtkScalarBarActor> scalarBar;
    vtkNew<vtkTextActor> actor_text;
    vtkNew<vtkTextActor> actor_text_title;
    vtkNew<vtkActor> rectangleActor;
    vtkNew<vtkActor> actor_wind;

private:
    vtkNew<vtkLookupTable> hueLut_count, hueLut_J, hueLut_Jpinv, hueLut_Q, hueLut;

    // main grid
    vtkNew<vtkStructuredGrid> structuredGrid_main;
    vtkNew<vtkDataSetMapper> mapper_grid_main;
    vtkNew<vtkPoints> grid_main_points;
    vtkNew<vtkUnsignedCharArray> pts_colors;

    // frame around the grid
    vtkNew<vtkPolyLine> polyline;
    vtkNew<vtkPolyData> rectanglePolyData;
    vtkNew<vtkPoints> rectanglePoints;
    vtkNew<vtkCellArray> rectangleLines;
    vtkNew<vtkPolyDataMapper> rectangleMapper;

    // lat/lon grid
    vtkNew<vtkStructuredGrid> structuredGrid_lat_lon;
    vtkNew<vtkDataSetMapper> mapper_grid_lat_lon;
    vtkNew<vtkPoints> grid_lat_lon_points;
    vtkNew<vtkFloatArray> visualized_values_grid;

    void populateLut(const float lutArray[][3], const int size, vtkNew<vtkLookupTable> &table);
    void interpolateLut(const float lutArray[][3], const int size, vtkNew<vtkLookupTable> &table);

    // wind visualization
    vtkNew<vtkPoints> points_wind_vector;
    vtkNew<vtkDoubleArray> vectors_values;
    vtkNew<vtkPolyData> polyData_wind;
    vtkNew<vtkArrowSource> arrowSource;
    vtkNew<vtkGlyph2D> glyphFilter;
    vtkNew<vtkPolyDataMapper> vectorMapper_wind;

    static constexpr int numwRows = 4;
    static constexpr int numwCols = 5;

static constexpr float lutSpecialJ[5][3] = {
        {0.342992, 0.650614, 0.772702},//10
        {0.688376, 0.931066, 0.963615},//7

        {0.836049, 0.882901, 0.85903},//5

        {0xee/255.,0x6e/255.,0xba/255.},

        {0x94/255.,0x00/255.,0x58/255.}
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


static constexpr float naturalRidges[7][3] = {
    {0x03/255.,0x36/255.,0xb3/255.},
    {0x27/255.,0x5d/255.,0x96/255.},    //3
    {0x80/255.,0xac/255.,0xd9/255.},    //1

    {0x77/255.,0x9a/255.,0xae/255.},  // crushed 9

    {0xdf/255.,0xa7/255.,0xac/255.},             // -1
    {0xc2/255.,0x66/255.,0x6e/255.},             // -2
    {0xb7/255.,0x24/255.,0x30/255.}

};

static constexpr float lutSpecialQ[9][3] = {
    {1., 0.997467, 0.914244},
    {0.89379, 0.887433, 0.810234},
    {0.788048, 0.780999, 0.707756},
    {0.685244, 0.687897, 0.609009},
    {0.58726, 0.614242, 0.515383},
    {0.495398, 0.562531, 0.427457},
    {0.416896, 0.521963, 0.367477},
    {0.360777, 0.485681, 0.367074},
    {0.323764, 0.451244, 0.413069}
};

static Eigen::Vector3f interpolateColor(const float colorArray[][3], int nColors, float value);

};
#endif
