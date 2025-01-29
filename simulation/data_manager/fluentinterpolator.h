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
    vtkNew<vtkDataSetMapper> mapper;
    vtkNew<vtkUnstructuredGrid> grid;
    vtkNew<vtkLookupTable> lut;


private:
    std::string geometryFilePrefix;
//    vtkNew<vtkFLUENTCFFReader> fluentReader;
    vtkNew<vtkFLUENTCFFCustomReader> fluentReader;
};

#endif // FLUENTINTERPOLATOR_H
