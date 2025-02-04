#ifndef FLUENTINTERPOLATOR_H
#define FLUENTINTERPOLATOR_H

#include <vtkNew.h>
//#include <vtkFLUENTCFFReader.h>
#include "vtkFLUENTCFFCustomReader.h"

#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkLookupTable.h>
#include <vtkCellData.h>
#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkProbeFilter.h>
#include <vtkImageData.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkTransformFilter.h>
#include <vtkTransform.h>

#include <string>
#include <utility>

#include <Eigen/Core>

#include "parameters_sim.h"

class FluentInterpolator
{
public:
    FluentInterpolator();
    SimParams *prms;
    bool is_initialized = false;

    void ScanDirectory(std::string geometryFile);

    void TestLoad(double scale, double ox, double oy);

    bool SetTime(double t);

    Eigen::Vector2f getInterpolation(int i, int j);

    float *getData(int frame, int subIdx);


    int file_count, interval_size;   // files from the scanned directory
    int currentFrame = -1;

    vtkNew<vtkActor> actor;
    vtkNew<vtkActor> actor_original;

    vtkNew<vtkDataSetMapper> mapper;
    vtkNew<vtkUnstructuredGrid> grid;
    vtkNew<vtkLookupTable> lut;
    vtkNew<vtkDataSetMapper> probeMapper;

private:
    double currentTime, position;

    std::string geometryFilePrefix;
    vtkNew<vtkFLUENTCFFCustomReader> fluentReader;
    vtkNew<vtkTransform> transform;
    vtkNew<vtkTransformFilter> transformFilter;
    vtkNew<vtkCellDataToPointData> filter_cd2pd;
    vtkNew<vtkImageData> imageData;
    vtkNew<vtkProbeFilter> probeFilter;

    std::vector<float> _data;

    void LoadDataFrame(int frame, float *U, float *V);

};

#endif // FLUENTINTERPOLATOR_H
