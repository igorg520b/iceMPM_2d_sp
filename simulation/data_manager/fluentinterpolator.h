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

#include <string>

class FluentInterpolator
{
public:
    FluentInterpolator();

    void ScanDirectory(std::string geometryFile);
    void LoadDataFrame(int frame);

    void TestLoad();


    int file_count, interval;

    vtkNew<vtkActor> actor;
    vtkNew<vtkActor> actor_original;

    vtkNew<vtkDataSetMapper> mapper;
    vtkNew<vtkUnstructuredGrid> grid;
    vtkNew<vtkLookupTable> lut;

    vtkNew<vtkCellDataToPointData> filter_cd2pd;
    vtkNew<vtkImageData> imageData;
    vtkNew<vtkProbeFilter> probeFilter;

    vtkNew<vtkDataSetMapper> probeMapper;


private:
    std::string geometryFilePrefix;
//    vtkNew<vtkFLUENTCFFReader> fluentReader;
    vtkNew<vtkFLUENTCFFCustomReader> fluentReader;

    std::vector<float> flatData_U, flatData_V;

};

#endif // FLUENTINTERPOLATOR_H
