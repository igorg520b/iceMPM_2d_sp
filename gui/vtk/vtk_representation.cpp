#include "vtk_representation.h"
#include "model.h"
#include "parameters_sim.h"
#include "gpu_partition.h"
//#include <omp.h>
#include <algorithm>
#include <iostream>
#include <spdlog/spdlog.h>

icy::VisualRepresentation::VisualRepresentation()
{
    int nLut = sizeof lutArrayTemperatureAdj / sizeof lutArrayTemperatureAdj[0];
    hueLut_temperature->SetNumberOfTableValues(nLut);
    for ( int i=0; i<nLut; i++)
        hueLut_temperature->SetTableValue(i, lutArrayTemperatureAdj[i][0],
                              lutArrayTemperatureAdj[i][1],
                              lutArrayTemperatureAdj[i][2], 1.0);


    nLut = sizeof(lutSouthwest)/sizeof lutSouthwest[0];
    hueLut_Southwest->SetNumberOfTableValues(nLut);
    for ( int i=0; i<nLut; i++)
        hueLut_Southwest->SetTableValue(i, lutSouthwest[i][0],
                                          lutSouthwest[i][1],
                                          lutSouthwest[i][2], 1.0);


    nLut = sizeof lutArrayPastel / sizeof lutArrayPastel[0];
    hueLut_pastel->SetNumberOfTableValues(nLut);
    for ( int i=0; i<nLut; i++)
        hueLut_pastel->SetTableValue(i, lutArrayPastel[i][0],
                              lutArrayPastel[i][1],
                              lutArrayPastel[i][2], 1.0);
    hueLut_pastel->SetTableRange(0,nLut-1);

    nLut = sizeof lutArrayMPMColors / sizeof lutArrayMPMColors[0];
    lutMPM->SetNumberOfTableValues(nLut);
    for ( int i=0; i<nLut; i++)
        lutMPM->SetTableValue(i, lutArrayMPMColors[i][0],
                              lutArrayMPMColors[i][1],
                              lutArrayMPMColors[i][2], 1.0);

    hueLut_four->SetNumberOfColors(5);
    hueLut_four->SetTableValue(0, 0.3, 0.3, 0.3);
    hueLut_four->SetTableValue(1, 1.0, 0, 0);
    hueLut_four->SetTableValue(2, 0, 1.0, 0);
    hueLut_four->SetTableValue(3, 0.2, 0.2, 0.85);
    hueLut_four->SetTableValue(4, 0, 0.5, 0.5);
    hueLut_four->SetTableRange(0,4);


    indenterMapper->SetInputConnection(indenterSource->GetOutputPort());
    actor_indenter->SetMapper(indenterMapper);

    indenterSource->GeneratePolygonOff();
    indenterSource->SetNumberOfSides(50);

    indenterMapper->SetInputConnection(indenterSource->GetOutputPort());
    actor_indenter->SetMapper(indenterMapper);
    actor_indenter->GetProperty()->LightingOff();
    actor_indenter->GetProperty()->EdgeVisibilityOn();
    actor_indenter->GetProperty()->VertexVisibilityOff();
    actor_indenter->GetProperty()->SetColor(0.1,0.1,0.1);
    actor_indenter->GetProperty()->SetEdgeColor(90.0/255.0, 90.0/255.0, 97.0/255.0);
    actor_indenter->GetProperty()->ShadingOff();
    actor_indenter->GetProperty()->SetInterpolationToFlat();
    actor_indenter->PickableOff();
    actor_indenter->GetProperty()->SetLineWidth(3);


    // points
    points_polydata->SetPoints(points);
    points_polydata->GetPointData()->AddArray(visualized_values);
    visualized_values->SetName("visualized_values");
    points_polydata->GetPointData()->SetActiveScalars("visualized_values");

    points_filter->SetInputData(points_polydata);
    points_filter->Update();

    points_mapper->SetInputData(points_filter->GetOutput());
    points_mapper->UseLookupTableScalarRangeOn();
    points_mapper->SetLookupTable(lutMPM);

    actor_points->SetMapper(points_mapper);
    actor_points->GetProperty()->SetPointSize(2);
    actor_points->GetProperty()->SetVertexColor(1,0,0);
    actor_points->GetProperty()->SetColor(0,0,0);
    actor_points->GetProperty()->LightingOff();
    actor_points->GetProperty()->ShadingOff();
    actor_points->GetProperty()->SetInterpolationToFlat();
    actor_points->PickableOff();

    grid_mapper->SetInputData(structuredGrid);
//    grid_mapper->SetLookupTable(hueLut);

    actor_grid->SetMapper(grid_mapper);
    actor_grid->GetProperty()->SetEdgeVisibility(true);
    actor_grid->GetProperty()->SetEdgeColor(0.8,0.8,0.8);
    actor_grid->GetProperty()->LightingOff();
    actor_grid->GetProperty()->ShadingOff();
    actor_grid->GetProperty()->SetInterpolationToFlat();
    actor_grid->PickableOff();
    actor_grid->GetProperty()->SetColor(0.98,0.98,0.98);
//    actor_grid->GetProperty()->SetRepresentationToWireframe();

    // scalar bar
    scalarBar->SetLookupTable(lutMPM);
    scalarBar->SetMaximumWidthInPixels(130);
    scalarBar->SetBarRatio(0.07);
    scalarBar->SetMaximumHeightInPixels(200);
    scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();
    scalarBar->GetPositionCoordinate()->SetValue(0.01,0.015, 0.0);
    scalarBar->SetLabelFormat("%.1e");
    scalarBar->GetLabelTextProperty()->BoldOff();
    scalarBar->GetLabelTextProperty()->ItalicOff();
    scalarBar->GetLabelTextProperty()->ShadowOff();
    scalarBar->GetLabelTextProperty()->SetColor(0.1,0.1,0.1);

    // text
    vtkTextProperty* txtprop = actorText->GetTextProperty();
    txtprop->SetFontFamilyToArial();
    txtprop->BoldOff();
    txtprop->SetFontSize(14);
    txtprop->ShadowOff();
    txtprop->SetColor(0,0,0);
    actorText->SetDisplayPosition(500, 30);

}



void icy::VisualRepresentation::SynchronizeTopology()
{
    model->accessing_point_data.lock();

    spdlog::info("SynchronizeTopology()");

    points->SetNumberOfPoints(model->prms.nPtsTotal);
    visualized_values->SetNumberOfValues(model->prms.nPtsTotal);
    indenterSource->SetRadius(model->prms.IndDiameter/2.f);


    int gx = model->prms.GridXTotal;
    int gy = model->prms.GridY;
    float &h = model->prms.cellsize;
    structuredGrid->SetDimensions(gx, gy, 1);

    grid_points->SetNumberOfPoints(gx*gy);
    for(int idx_y=0; idx_y<gy; idx_y++)
        for(int idx_x=0; idx_x<gx; idx_x++)
        {
            float x = idx_x * h;
            float y = idx_y * h;
            double pt_pos[3] {x, y, -2.0};
            grid_points->SetPoint((vtkIdType)(idx_x+idx_y*gx), pt_pos);
        }
    structuredGrid->SetPoints(grid_points);


    model->accessing_point_data.unlock();
    SynchronizeValues();
}


void icy::VisualRepresentation::SynchronizeValues()
{
    model->accessing_point_data.lock();

    double sim_time = model->prms.SimulationTime;

    // spdlog::info("SynchronizeValues() npts {}", model->prms.nPtsTotal);
    double indenter_x = model->prms.indenter_x;
    double indenter_y = model->prms.indenter_y;
    indenterSource->SetCenter(indenter_x, indenter_y, 1);


    unsigned activePtsCount = 0;
    for(int i=0;i<model->gpu.hssoa.size;i++)
    {
        SOAIterator s = model->gpu.hssoa.begin()+i;
        if(s->getDisabledStatus()) continue;
        Eigen::Vector2f pos = s->getPos(model->prms.cellsize);
        points->SetPoint((vtkIdType)activePtsCount, pos[0], pos[1], 0);
//        spdlog::info("setting point {}; {} - {}", activePtsCount, pos[0], pos[1]);
        activePtsCount++;
    }
    if(activePtsCount != model->prms.nPtsTotal) throw std::runtime_error("SynchronizeValues() point count mismatch (pos)");
    points->Modified();
    actor_points->GetProperty()->SetPointSize(model->prms.ParticleViewSize);
    points_filter->Update();
    double range = std::pow(10,ranges[VisualizingVariable]);
    double centerVal = 0;



    if(VisualizingVariable == VisOpt::status)
    {
        scalarBar->VisibilityOn();
        points_mapper->ScalarVisibilityOn();
        points_mapper->SetColorModeToMapScalars();
        points_mapper->UseLookupTableScalarRangeOn();
        points_mapper->SetLookupTable(hueLut_four);
        scalarBar->SetLookupTable(hueLut_four);
        visualized_values->SetNumberOfValues(model->prms.nPtsTotal);
        activePtsCount = 0;
        for(int i=0;i<model->gpu.hssoa.size;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            if(s->getDisabledStatus()) continue;
            bool isCrushed = s->getCrushedStatus();
            float value = 0;
            if(isCrushed) value = 1;
            visualized_values->SetValue((vtkIdType)activePtsCount++, (float)value);
        }
        if(activePtsCount != model->prms.nPtsTotal) throw std::runtime_error("SynchronizeValues() point count mismatch");
        visualized_values->Modified();
    }
    else if(VisualizingVariable == VisOpt::Jp_inv)
    {
        scalarBar->VisibilityOn();
        points_mapper->ScalarVisibilityOn();
        points_mapper->SetColorModeToMapScalars();
        points_mapper->UseLookupTableScalarRangeOn();
        points_mapper->SetLookupTable(lutMPM);
        scalarBar->SetLookupTable(lutMPM);
        lutMPM->SetTableRange(centerVal-range, centerVal+range);
        visualized_values->SetNumberOfValues(model->prms.nPtsTotal);
        activePtsCount = 0;
        for(int i=0;i<model->gpu.hssoa.size;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            if(s->getDisabledStatus()) continue;
            double value = s->getValue(SimParams::idx_Jp_inv)-1;
            visualized_values->SetValue((vtkIdType)activePtsCount++, (float)value);
        }
        if(activePtsCount != model->prms.nPtsTotal) throw std::runtime_error("SynchronizeValues() point count mismatch");
        visualized_values->Modified();
    }
    else if(VisualizingVariable == VisOpt::grains)
    {
        scalarBar->VisibilityOn();
        points_mapper->ScalarVisibilityOn();
        points_mapper->SetColorModeToMapScalars();
        points_mapper->UseLookupTableScalarRangeOn();
        points_mapper->SetLookupTable(hueLut_pastel);
        scalarBar->SetLookupTable(hueLut_pastel);

        visualized_values->SetNumberOfValues(model->prms.nPtsTotal);
        activePtsCount = 0;
        for(int i=0;i<model->gpu.hssoa.size;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            if(s->getDisabledStatus()) continue;
            uint16_t grain = s->getGrain()%40;
            bool isCrushed = s->getCrushedStatus();
            if(isCrushed) grain = 41;
            visualized_values->SetValue((vtkIdType)activePtsCount++, (float)grain);
        }
        if(activePtsCount != model->prms.nPtsTotal) throw std::runtime_error("SynchronizeValues() point count mismatch");
        visualized_values->Modified();
    }
    else if(VisualizingVariable == VisOpt::velocity)
    {
        scalarBar->VisibilityOn();
        points_mapper->ScalarVisibilityOn();
        points_mapper->SetColorModeToMapScalars();
        points_mapper->UseLookupTableScalarRangeOn();
        points_mapper->SetLookupTable(hueLut_Southwest);
        scalarBar->SetLookupTable(hueLut_Southwest);
//        hueLut_temperature->SetTableRange(0, 7.0);
        hueLut_Southwest->SetRange(0, 2.0);

        //points_mapper->SetScalarRange(11,-1);
        visualized_values->SetNumberOfValues(model->prms.nPtsTotal);
        activePtsCount = 0;
        for(int i=0;i<model->gpu.hssoa.size;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            if(s->getDisabledStatus()) continue;
            Eigen::Vector2f vel = s->getVelocity();
            float value = (float)vel.norm();
            if(value>=1.9) value=1.9;
            visualized_values->SetValue((vtkIdType)activePtsCount++, (float)value);
        }
        if(activePtsCount != model->prms.nPtsTotal) throw std::runtime_error("SynchronizeValues() point count mismatch");
        visualized_values->Modified();
    }
    else
    {
        points_mapper->ScalarVisibilityOff();
        scalarBar->VisibilityOff();
    }
    model->accessing_point_data.unlock();

}


void icy::VisualRepresentation::ChangeVisualizationOption(int option)
{
    VisualizingVariable = (VisOpt)option;
    SynchronizeTopology();
}

