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

#include "parameters_sim.h"

class FluentInterpolator
{
public:
    FluentInterpolator();
    SimParams *prms;

    void ScanDirectory(std::string geometryFile);

    void TestLoad(double scale, double ox, double oy);

    bool SetTime(double t);



    int file_count, interval_size;   // files from the scanned directory
    std::vector<float> flatData1_U, flatData1_V;  // accessible by the model to transfer to GPU
    std::vector<float> flatData2_U, flatData2_V;  // value of interest is interpolated between (1) and (2)



    vtkNew<vtkActor> actor;
    vtkNew<vtkActor> actor_original;

    vtkNew<vtkDataSetMapper> mapper;
    vtkNew<vtkUnstructuredGrid> grid;
    vtkNew<vtkLookupTable> lut;

    vtkNew<vtkCellDataToPointData> filter_cd2pd;
    vtkNew<vtkImageData> imageData;
    vtkNew<vtkProbeFilter> probeFilter;

    vtkNew<vtkDataSetMapper> probeMapper;
    vtkNew<vtkTransformFilter> transformFilter;
    vtkNew<vtkTransform> transform;


private:
    int gx, gy; // modelled area grid resolution
    int ox, oy; // offset of modelled area from visual grid
    int vgx, vgy; // visual grid dimensions
    double scale, offsetX, offsetY;
    int currentFrame = -1;
    double currentTime, position;

    std::string geometryFilePrefix;
    vtkNew<vtkFLUENTCFFCustomReader> fluentReader;

    void LoadDataFrame(int frame, std::vector<float> &flatData_U, std::vector<float> &flatData_V);

};

#endif // FLUENTINTERPOLATOR_H
