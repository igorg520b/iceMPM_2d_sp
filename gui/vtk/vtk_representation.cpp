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


    pts_colors->SetNumberOfComponents(3);
    pts_colors->SetName("pts_colors");
    points_polydata->GetPointData()->AddArray(pts_colors);

    // GRID
    grid_mapper->SetInputData(structuredGrid);
//    grid_mapper->SetLookupTable(hueLut);

    actor_grid->SetMapper(grid_mapper);
    actor_grid->GetProperty()->SetEdgeVisibility(true);
    actor_grid->GetProperty()->SetEdgeColor(0.2,0.2,0.2);
    actor_grid->GetProperty()->LightingOff();
    actor_grid->GetProperty()->ShadingOff();
    actor_grid->GetProperty()->SetInterpolationToFlat();
    actor_grid->PickableOff();
    actor_grid->GetProperty()->SetColor(0.1,0.1,0.1);
    actor_grid->GetProperty()->SetRepresentationToWireframe();

    // grid2

    grid_colors->SetName("grid_colors");
    grid_colors->SetNumberOfComponents(3);

    structuredGrid2->GetCellData()->AddArray(grid_colors);
    structuredGrid2->GetCellData()->SetActiveScalars("grid_colors");

    grid_mapper2->SetInputData(structuredGrid2);
//    grid_mapper2->SetLookupTable(hueLut_four);
//    grid_mapper2->UseLookupTableScalarRangeOn();
    grid_mapper2->SetColorModeToDirectScalars();


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


    // uniform grid
    mapper_uniformgrid->SetInputData(uniformGrid);
    actor_uniformgrid->SetMapper(mapper_uniformgrid);
}



void icy::VisualRepresentation::SynchronizeTopology()
{
    model->accessing_point_data.lock();

    spdlog::info("SynchronizeTopology()");

    const int nPts = model->gpu.hssoa.size;
    points->SetNumberOfPoints(nPts);
    visualized_values->SetNumberOfValues(nPts);
    pts_colors->SetNumberOfValues(nPts*3);

    int gx = model->prms.GridXTotal;
    int gy = model->prms.GridYTotal;

    const int imgx = model->prms.InitializationImageSizeX;
    const int imgy = model->prms.InitializationImageSizeY;
    const double h = model->prms.cellsize;

    uniformGrid->SetDimensions(imgx, imgy, 1);
    uniformGrid->SetSpacing(h, h, 1.0);
    uniformGrid->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    for(int i=0;i<imgx;i++)
        for(int j=0;j<imgy;j++)
        {
            unsigned char* pixel = static_cast<unsigned char*>(uniformGrid->GetScalarPointer(i, j, 0));
            int idx = (i+imgx*j)*3;
            pixel[0] = model->gpu.original_image_colors_rgb[idx+0];
            pixel[1] = model->gpu.original_image_colors_rgb[idx+1];
            pixel[2] = model->gpu.original_image_colors_rgb[idx+2];
        }
//    uniformGrid->Modified();
    mapper_uniformgrid->Update();

/*
    // GRID1 now shows the parallels and meridians
    int nGridLat = model->wind_interpolator.extentLat;
    int nGridLon = model->wind_interpolator.extentLon;

    double lami = model->wind_interpolator.LatMin;
    double lamx = model->wind_interpolator.LatMax;
    double lomi = model->wind_interpolator.LonMin;
    double lomx = model->wind_interpolator.LonMax;
    auto latToGrid = [&](float lat) { return h*gy*(lat-lami)/(lamx-lami); };
    auto lonToGrid = [&](float lon) { return h*gx*(lon-lomi)/(lomx-lomi); };

    if(model->wind_interpolator.isInitialized)
    {
        structuredGrid->SetDimensions(nGridLon, nGridLat, 1);
        grid_points->SetNumberOfPoints(nGridLon*nGridLat);
        for(int idx_y=0; idx_y<nGridLat; idx_y++)
            for(int idx_x=0; idx_x<nGridLon; idx_x++)
            {
                float lon = model->wind_interpolator.gridLonMin + idx_x*WindInterpolator::gridCellSize;
                float lat = model->wind_interpolator.gridLatMin + idx_y*WindInterpolator::gridCellSize;
                double x = lonToGrid(lon);
                double y = latToGrid(lat);
                double pt_pos[3] {x, y, -0.5};
                grid_points->SetPoint((vtkIdType)(idx_x+idx_y*nGridLon), pt_pos);
            }
        structuredGrid->SetPoints(grid_points);
    }
  */


    /*

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

//    visualized_values_grid->SetNumberOfValues(gx*gy);
    grid_colors->SetNumberOfValues(gx*gy*3);
    std::vector<uint8_t> &gb = model->gpu.hssoa.grid_status_buffer;
    if(gb.size() > 0)
    {
        spdlog::info("gb size {}; gx {}; gy {}",gb.size(), gx, gy);
        for(int idx_y=0; idx_y<gy; idx_y++)
            for(int idx_x=0; idx_x<gx; idx_x++)
            {
                int idx_visual = idx_x + idx_y*gx;
                int idx_storage = idx_y + idx_x*gy;

                //int value = gb[idx_storage]+5;
                //visualized_values_grid->SetValue(idx_visual, value);
                uint8_t r = model->gpu.hssoa.grid_colors_rgb[idx_storage*3+0];
                uint8_t g = model->gpu.hssoa.grid_colors_rgb[idx_storage*3+1];
                uint8_t b = model->gpu.hssoa.grid_colors_rgb[idx_storage*3+2];

                uint8_t water = model->gpu.hssoa.grid_status_buffer[idx_storage];
                if(water) {r = 29; g=41; b=58;}

                grid_colors->SetValue(idx_visual*3+0, r);
                grid_colors->SetValue(idx_visual*3+1, g);
                grid_colors->SetValue(idx_visual*3+2, b);
            }

        grid_colors->Modified();
    }
*/


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
        points->SetPoint((vtkIdType)i, pos[0], pos[1], +1.0);
    }
    points->Modified();
    actor_points->GetProperty()->SetPointSize(model->prms.ParticleViewSize);
    double range = std::pow(10,ranges[VisualizingVariable]);
    double centerVal = 0;

    actor_points->VisibilityOn();


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
            float value = 1;
            SOAIterator s = model->gpu.hssoa.begin()+i;
            if(s->getDisabledStatus()) value = 3;
            else if(s->getCrushedStatus()) value = 0;
            else if(s->getWeakenedStatus()) value = 2;
            visualized_values->SetValue((vtkIdType)i, (float)value);
        }
        points_polydata->GetPointData()->SetActiveScalars("visualized_values");
        visualized_values->Modified();
        grid_mapper2->SetLookupTable(hueLut_four);
    }

    else if(VisualizingVariable == VisOpt::wind_u)
    {
        /*
        int gx = model->prms.GridXTotal;
        int gy = model->prms.GridY;
        actor_points->VisibilityOff();

        model->wind_interpolator.setTime(wind_visualization_time);
        float tb = model->wind_interpolator.interpolationCoeffFromTime(wind_visualization_time);

        for(int idx_y=0; idx_y<gy; idx_y++)
            for(int idx_x=0; idx_x<gx; idx_x++)
            {
                int idx_visual = idx_x + idx_y*gx;

                float lat = model->prms.LatMin + (model->prms.LatMax-model->prms.LatMin)*(float)idx_y/(float)gy;
                float lon = model->prms.LonMin + (model->prms.LonMax-model->prms.LonMin)*(float)idx_x/(float)gx;

                Eigen::Vector2f wv = model->wind_interpolator.interpolationResult(lat, lon, tb);
                float value = wv.x();
                visualized_values_grid->SetValue(idx_visual, value);
            }

        grid_mapper2->SetLookupTable(hueLut_pressure);
        hueLut_pressure->SetTableRange(-range, range);

        visualized_values_grid->Modified();
        scalarBar->VisibilityOn();
        scalarBar->SetLookupTable(hueLut_pressure);
*/
    }
    else if(VisualizingVariable == VisOpt::wind_v)
    {
        /*
        int gx = model->prms.GridXTotal;
        int gy = model->prms.GridY;
        actor_points->VisibilityOff();
        model->wind_interpolator.setTime(wind_visualization_time);
        float tb = model->wind_interpolator.interpolationCoeffFromTime(wind_visualization_time);

        for(int idx_y=0; idx_y<gy; idx_y++)
            for(int idx_x=0; idx_x<gx; idx_x++)
            {
                int idx_visual = idx_x + idx_y*gx;

                float lat = model->prms.LatMin + (model->prms.LatMax-model->prms.LatMin)*(float)idx_y/(float)gy;
                float lon = model->prms.LonMin + (model->prms.LonMax-model->prms.LonMin)*(float)idx_x/(float)gx;

                Eigen::Vector2f wv = model->wind_interpolator.interpolationResult(lat, lon, tb);
                float value = wv.y();
                visualized_values_grid->SetValue(idx_visual, value);
            }

        grid_mapper2->SetLookupTable(hueLut_pressure);
        hueLut_pressure->SetTableRange(-range, range);

        visualized_values_grid->Modified();
        scalarBar->VisibilityOn();
        scalarBar->SetLookupTable(hueLut_pressure);
*/
    }
    else if(VisualizingVariable == VisOpt::wind_norm)
    {
        /*
        int gx = model->prms.GridXTotal;
        int gy = model->prms.GridY;
        actor_points->VisibilityOff();

        model->wind_interpolator.setTime(wind_visualization_time);
        float tb = model->wind_interpolator.interpolationCoeffFromTime(wind_visualization_time);

        for(int idx_y=0; idx_y<gy; idx_y++)
            for(int idx_x=0; idx_x<gx; idx_x++)
            {
                int idx_visual = idx_x + idx_y*gx;

                float lat = model->prms.LatMin + (model->prms.LatMax-model->prms.LatMin)*(float)idx_y/(float)gy;
                float lon = model->prms.LonMin + (model->prms.LonMax-model->prms.LonMin)*(float)idx_x/(float)gx;

                Eigen::Vector2f wv = model->wind_interpolator.interpolationResult(lat, lon, tb);
                float value = wv.norm();
                visualized_values_grid->SetValue(idx_visual, value);
            }

        grid_mapper2->SetLookupTable(hueLut_Southwest);
        hueLut_Southwest->SetTableRange(0, range);

        scalarBar->VisibilityOn();
        scalarBar->SetLookupTable(hueLut_Southwest);

        visualized_values_grid->Modified();
*/
    }
    else if(VisualizingVariable == VisOpt::color)
    {
        scalarBar->VisibilityOff();
        points_mapper->ScalarVisibilityOn();

        points_mapper->SetColorModeToDirectScalars();
//        points_mapper->UseLookupTableScalarRangeOff();
        points_mapper->Modified();

        for(int i=0;i<nPts;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            int pt_idx = s->getValueInt(SimParams::integer_point_idx);
            uint32_t rgb = model->gpu.point_colors_rgb[pt_idx];
            uint8_t r = (rgb >> 16) & 0xff;
            uint8_t g = (rgb >> 8) & 0xff;
            uint8_t b = rgb & 0xff;
            pts_colors->SetValue((vtkIdType)(i*3+0), r);
            pts_colors->SetValue((vtkIdType)(i*3+1), g);
            pts_colors->SetValue((vtkIdType)(i*3+2), b);
        }
        points_polydata->GetPointData()->SetActiveScalars("pts_colors");
        pts_colors->Modified();

        points_polydata->Modified();
        points_filter->Update();
        grid_mapper2->SetLookupTable(hueLut_four);
    }
    else if(VisualizingVariable == VisOpt::special)
    {
        scalarBar->VisibilityOff();
        points_mapper->ScalarVisibilityOn();

        points_mapper->SetColorModeToDirectScalars();
        points_mapper->Modified();

        for(int i=0;i<nPts;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            int pt_idx = s->getValueInt(SimParams::integer_point_idx);
            uint32_t rgb = model->gpu.point_colors_rgb[pt_idx];
            uint8_t _r = (rgb >> 16) & 0xff;
            uint8_t _g = (rgb >> 8) & 0xff;
            uint8_t _b = rgb & 0xff;

            float r = _r/255.;
            float g = _g/255.;
            float b = _b/255.;

            if(s->getCrushedStatus())
            {
                double value = s->getValue(SimParams::idx_Jp_inv)-1;
                interpolateColor(naturalRidges, 7, (value)*5+0.5, r, g, b);
            }

            pts_colors->SetValue((vtkIdType)(i*3+0), (unsigned char)(r*255));
            pts_colors->SetValue((vtkIdType)(i*3+1), (unsigned char)(g*255));
            pts_colors->SetValue((vtkIdType)(i*3+2), (unsigned char)(b*255));
        }
        points_polydata->GetPointData()->SetActiveScalars("pts_colors");
        pts_colors->Modified();

        points_polydata->Modified();
        points_filter->Update();
        grid_mapper2->SetLookupTable(hueLut_four);
    }
    else if(VisualizingVariable == VisOpt::Jp_inv)
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
            double value = s->getValue(SimParams::idx_Jp_inv)-1;
//            if(s->getDisabledStatus()) value = 0;
//            else if(!s->getCrushedStatus()) value = -10;
            visualized_values->SetValue((vtkIdType)i, (float)value);
        }
        points_polydata->GetPointData()->SetActiveScalars("visualized_values");
        visualized_values->Modified();
        grid_mapper2->SetLookupTable(hueLut_four);
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
            // else if(s->getCrushedStatus()) value = 2*range;
            visualized_values->SetValue((vtkIdType)i, value);
        }
        points_polydata->GetPointData()->SetActiveScalars("visualized_values");
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
            //else if(s->getCrushedStatus()) value = 2*range;
            visualized_values->SetValue((vtkIdType)i, value);
        }
        points_polydata->GetPointData()->SetActiveScalars("visualized_values");
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
        points_polydata->GetPointData()->SetActiveScalars("visualized_values");
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
        points_polydata->GetPointData()->SetActiveScalars("visualized_values");
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
        points_polydata->GetPointData()->SetActiveScalars("visualized_values");
        visualized_values->Modified();
    }
    else if(VisualizingVariable == VisOpt::boundary)
    {
        actor_points->VisibilityOff();
    }
    else if(VisualizingVariable == VisOpt::strength)
    {
        scalarBar->VisibilityOn();
        points_mapper->ScalarVisibilityOn();
        points_mapper->SetColorModeToMapScalars();
        points_mapper->UseLookupTableScalarRangeOn();
        points_mapper->SetLookupTable(hueLut_pressure);
        scalarBar->SetLookupTable(hueLut_pressure);
        hueLut_pressure->SetTableRange(0.5, 1);
        for(int i=0;i<nPts;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            float value = s->getValue(SimParams::idx_initial_strength);
            //else if(s->getCrushedStatus()) value = 2*range;
            visualized_values->SetValue((vtkIdType)i, value);
        }
        points_polydata->GetPointData()->SetActiveScalars("visualized_values");
        visualized_values->Modified();
    }
    else
    {
        points_mapper->ScalarVisibilityOff();
        scalarBar->VisibilityOff();
    }
    points_filter->Update();

    model->accessing_point_data.unlock();

}


void icy::VisualRepresentation::ChangeVisualizationOption(int option)
{
    VisualizingVariable = (VisOpt)option;
    SynchronizeTopology();
}

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


void icy::VisualRepresentation::interpolateColor(const float colorArray[][3],
                                                 int nColors, float value, float& r_out, float& g_out, float& b_out) {
    if (nColors < 2) {
        throw std::invalid_argument("Color array must have at least two colors for interpolation.");
    }

    // Clamp value to range [0, 1]
    value = std::clamp(value, 0.0f, 1.0f);

    // Find the segment
    float segmentSize = 1.0f / (nColors - 1); // Divide [0, 1] into nColors-1 segments
    int segmentIndex = static_cast<int>(value / segmentSize); // Determine which segment we're in
    segmentIndex = std::min(segmentIndex, nColors - 2); // Ensure index doesn't go out of bounds

    // Compute local position within the segment
    float localT = (value - segmentIndex * segmentSize) / segmentSize;

    // Interpolate between the two colors
    const float* color1 = colorArray[segmentIndex];
    const float* color2 = colorArray[segmentIndex + 1];
    r_out = (1 - localT) * color1[0] + localT * color2[0];
    g_out = (1 - localT) * color1[1] + localT * color2[1];
    b_out = (1 - localT) * color1[2] + localT * color2[2];
}
