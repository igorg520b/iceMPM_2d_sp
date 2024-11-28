#include "vtk_representation.h"
#include "model.h"
#include "parameters_sim.h"
#include "gpu_partition.h"
//#include <omp.h>
#include <algorithm>
#include <iostream>
#include <spdlog/spdlog.h>


void icy::VisualRepresentation::populateLut(const float lutArray[][3], const int size, vtkNew<vtkLookupTable> &table)
{
    table->SetNumberOfTableValues(size);
    for(int i=0;i<size;i++)
        table->SetTableValue(i,lutArray[i][0],lutArray[i][1],lutArray[i][2]);
    table->SetTableRange(0, size-1);

    table->SetRampToLinear();
}


void icy::VisualRepresentation::interpolateLut(const float lutArray[][3], const int size, vtkNew<vtkLookupTable> &table)
{
    const int m = 256;
    table->SetNumberOfTableValues(m);

    for (int i = 0; i < m; ++i) {
        float t = static_cast<float>(i) / (float)(m-1); // Normalize index to [0, 1]

        // Map t to the interval [0, N-1] for interpolation
        float scaledT = t * (size - 1);
        int lowerIdx = static_cast<int>(std::floor(scaledT)); // Lower bound
        int upperIdx = static_cast<int>(std::ceil(scaledT));  // Upper bound
        float localT = scaledT - lowerIdx; // Fractional position between indices

        // Linearly interpolate between colors at lowerIdx and upperIdx
        float r = (1.0f - localT) * lutArray[lowerIdx][0] + localT * lutArray[upperIdx][0];
        float g = (1.0f - localT) * lutArray[lowerIdx][1] + localT * lutArray[upperIdx][1];
        float b = (1.0f - localT) * lutArray[lowerIdx][2] + localT * lutArray[upperIdx][2];

        // Set the interpolated color in the lookup table
        table->SetTableValue(i, r, g, b, 1.0); // Alpha is set to 1.0 (opaque)
    }
}

icy::VisualRepresentation::VisualRepresentation()
{
    populateLut(lutArrayTemperatureAdj, (sizeof(lutArrayTemperatureAdj)/(sizeof(lutArrayTemperatureAdj[0]))), hueLut_temperature);
    populateLut(lutSouthwest, (sizeof(lutSouthwest)/(sizeof(lutSouthwest[0]))), hueLut_Southwest);
    populateLut(lutArrayPastel, (sizeof(lutArrayPastel)/(sizeof(lutArrayPastel[0]))), hueLut_pastel);
    populateLut(lutArrayMPMColors, (sizeof(lutArrayMPMColors)/(sizeof(lutArrayMPMColors[0]))), lutMPM);
    populateLut(lutSpecialSeven, (sizeof(lutSpecialSeven)/(sizeof(lutSpecialSeven[0]))), hueLut_four);

    interpolateLut(lutSpecialP, (sizeof(lutSpecialP)/(sizeof(lutSpecialP[0]))), hueLut_pressure);

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
//    actor_points->GetProperty()->SetColor(41./255,128./255,215./255);
    actor_points->GetProperty()->SetColor(151./255,188./255,215./255);
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

    // grid2
    visualized_values_grid->SetName("visualized_values_grid");
    structuredGrid2->GetCellData()->AddArray(visualized_values_grid);
    structuredGrid2->GetCellData()->SetActiveScalars("visualized_values_grid");

    grid_mapper2->SetInputData(structuredGrid2);
    grid_mapper2->SetLookupTable(hueLut_four);
    grid_mapper2->UseLookupTableScalarRangeOn();
    points_mapper->UseLookupTableScalarRangeOn();


    actor_grid2->SetMapper(grid_mapper2);
    actor_grid2->GetProperty()->SetEdgeVisibility(false);
    actor_grid2->GetProperty()->LightingOff();
    actor_grid2->GetProperty()->ShadingOff();
    actor_grid2->GetProperty()->SetInterpolationToFlat();
    actor_grid2->PickableOff();
    actor_grid2->GetProperty()->SetColor(0.98,0.98,0.98);

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

    const int nPts = model->gpu.hssoa.size;
    points->SetNumberOfPoints(nPts);
    visualized_values->SetNumberOfValues(nPts);

    int gx = model->prms.GridXTotal;
    int gy = model->prms.GridY;
    t_PointReal &h = model->prms.cellsize;
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

    // grid2
    int gx1 = gx+1;
    int gy1 = gy+1;
    structuredGrid2->SetDimensions(gx1, gy1, 1);
    grid_points2->SetNumberOfPoints(gx1*gy1);
    for(int idx_y=0; idx_y<=gy; idx_y++)
        for(int idx_x=0; idx_x<=gx; idx_x++)
        {
            float x = ((float)idx_x - 0.5) * h;
            float y = ((float)idx_y - 0.5) * h;
            double pt_pos[3] {x, y, -3.0};
            grid_points2->SetPoint((vtkIdType)(idx_x+idx_y*gx1), pt_pos);
        }
    structuredGrid2->SetPoints(grid_points2);

    visualized_values_grid->SetNumberOfValues(gx*gy);
    std::vector<uint8_t> &gb = model->gpu.hssoa.grid_status_buffer;
    if(gb.size() > 0)
    {
        spdlog::info("gb size {}; gx {}; gy {}",gb.size(), gx, gy);
        for(int idx_y=0; idx_y<gy; idx_y++)
            for(int idx_x=0; idx_x<gx; idx_x++)
            {
                int idx_visual = idx_x + idx_y*gx;
                int idx_storage = idx_y + idx_x*gy;
                int value = gb[idx_storage]+5;
                visualized_values_grid->SetValue(idx_visual, value);
            }

        visualized_values_grid->Modified();
    }



    model->accessing_point_data.unlock();
    SynchronizeValues();
}


void icy::VisualRepresentation::SynchronizeValues()
{
    model->accessing_point_data.lock();

    double sim_time = model->prms.SimulationTime;

    const int nPts = model->gpu.hssoa.size;
    for(int i=0;i<nPts;i++)
    {
        SOAIterator s = model->gpu.hssoa.begin()+i;
        PointVector2r pos = s->getPos(model->prms.cellsize);
        points->SetPoint((vtkIdType)i, pos[0], pos[1], 0);
    }
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
        for(int i=0;i<nPts;i++)
        {
            float value = 0;
            SOAIterator s = model->gpu.hssoa.begin()+i;
            if(s->getDisabledStatus()) value = 2;
            else if(s->getCrushedStatus()) value = 1;
            visualized_values->SetValue((vtkIdType)i, (float)value);
        }
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
        for(int i=0;i<nPts;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            double value = s->getValue(SimParams::idx_Jp_inv)-1;
            if(s->getDisabledStatus()) value = 0;
            visualized_values->SetValue((vtkIdType)i, (float)value);
        }
        visualized_values->Modified();
    }
    else if(VisualizingVariable == VisOpt::P)
    {
        scalarBar->VisibilityOn();
        points_mapper->ScalarVisibilityOn();
        points_mapper->SetColorModeToMapScalars();
        points_mapper->UseLookupTableScalarRangeOn();
        points_mapper->SetLookupTable(hueLut_pressure);
        scalarBar->SetLookupTable(hueLut_pressure);
        hueLut_pressure->SetTableRange(centerVal-range, centerVal+range);
        for(int i=0;i<nPts;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            float value = s->getValue(SimParams::idx_P);
            if(s->getDisabledStatus()) value = -2*range;
            else if(s->getCrushedStatus()) value = 2*range;
            visualized_values->SetValue((vtkIdType)i, value);
        }
        visualized_values->Modified();
    }
    else if(VisualizingVariable == VisOpt::Q)
    {
        scalarBar->VisibilityOn();
        points_mapper->ScalarVisibilityOn();
        points_mapper->SetColorModeToMapScalars();
        points_mapper->UseLookupTableScalarRangeOn();
        points_mapper->SetLookupTable(hueLut_pressure);
        scalarBar->SetLookupTable(hueLut_pressure);
        hueLut_pressure->SetTableRange(centerVal-range, centerVal+range);
        for(int i=0;i<nPts;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            float value = s->getValue(SimParams::idx_Q);
            if(s->getDisabledStatus()) value = -2*range;
            else if(s->getCrushedStatus()) value = 2*range;
            visualized_values->SetValue((vtkIdType)i, value);
        }
        visualized_values->Modified();
    }
    else if(VisualizingVariable == VisOpt::qp)
    {
        scalarBar->VisibilityOn();
        points_mapper->ScalarVisibilityOn();
        points_mapper->SetColorModeToMapScalars();
        points_mapper->UseLookupTableScalarRangeOn();
        points_mapper->SetLookupTable(hueLut_Southwest);
        scalarBar->SetLookupTable(hueLut_Southwest);
        hueLut_Southwest->SetTableRange(centerVal-range, centerVal+range);
        for(int i=0;i<nPts;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            float value = s->getValue(SimParams::idx_Qp)-1;
            if(s->getDisabledStatus()) value = -2*range;
            else if(s->getCrushedStatus()) value = 2*range;
            visualized_values->SetValue((vtkIdType)i, value);
        }
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
        for(int i=0;i<nPts;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            uint16_t grain = s->getGrain()%40;
            bool isCrushed = s->getCrushedStatus();
            if(s->getDisabledStatus()) grain = 41;
            if(isCrushed) grain = 41;
            visualized_values->SetValue((vtkIdType)i, (float)grain);
        }
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
        for(int i=0;i<nPts;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            PointVector2r vel = s->getVelocity();
            float value = (float)vel.norm();
            if(s->getDisabledStatus()) value = 1.9;
            if(value>=1.9) value=1.9;
            visualized_values->SetValue((vtkIdType)i, (float)value);
        }
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

