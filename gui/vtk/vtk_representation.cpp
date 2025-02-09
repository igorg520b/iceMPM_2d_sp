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
    // points
    pts_colors->SetNumberOfComponents(3);
    pts_colors->SetName("pts_colors");

    points_polydata->SetPoints(points);
    points_polydata->GetPointData()->AddArray(pts_colors);

    points_filter->SetInputData(points_polydata);
    points_filter->Update();

    points_mapper->SetInputData(points_filter->GetOutput());

    actor_points->SetMapper(points_mapper);
    actor_points->GetProperty()->SetPointSize(2);
    actor_points->GetProperty()->SetVertexColor(1,0,0);
    actor_points->GetProperty()->SetColor(151./255,188./255,215./255);
    actor_points->GetProperty()->LightingOff();
    actor_points->GetProperty()->ShadingOff();
    actor_points->GetProperty()->SetInterpolationToFlat();
    actor_points->PickableOff();


    // GRID
    grid_mapper->SetInputData(structuredGrid);

    actor_grid->SetMapper(grid_mapper);
    actor_grid->GetProperty()->SetEdgeVisibility(true);
    actor_grid->GetProperty()->SetEdgeColor(0.2,0.2,0.2);
    actor_grid->GetProperty()->LightingOff();
    actor_grid->GetProperty()->ShadingOff();
    actor_grid->GetProperty()->SetInterpolationToFlat();
    actor_grid->PickableOff();
    actor_grid->GetProperty()->SetColor(0.1,0.1,0.1);
    actor_grid->GetProperty()->SetRepresentationToWireframe();

    // scalar bar
    //scalarBar->SetLookupTable(lutMPM);
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
    spdlog::info("SynchronizeTopology()");

    const int nPts = model->gpu.hssoa.size;
    if(!nPts) return;

    model->accessing_point_data.lock();
    points->SetNumberOfPoints(nPts);
//    visualized_values->SetNumberOfValues(nPts);
    pts_colors->SetNumberOfValues(nPts*3);

    int ox = model->prms.ModeledRegionOffsetX;
    int oy = model->prms.ModeledRegionOffsetY;

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

    for(int i=0;i<gx;i++)
        for(int j=0;j<gy;j++)
        {
            int idx2 = (j + gy*i);
            if(model->gpu.grid_status_buffer[idx2])
            {
                // water color
                unsigned char* pixel = static_cast<unsigned char*>(uniformGrid->GetScalarPointer(i+ox, j+oy, 0));
                pixel[0] = 0x15;
                pixel[1] = 0x1f;
                pixel[2] = 0x2f;
            }
        }

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




    model->accessing_point_data.unlock();
    SynchronizeValues();
}


void icy::VisualRepresentation::SynchronizeValues()
{
    //spdlog::info("SynchronizeValues()");
    const int &_ox = model->prms.ModeledRegionOffsetX;
    const int &_oy = model->prms.ModeledRegionOffsetY;

    const int &gx = model->prms.GridXTotal;
    const int &gy = model->prms.GridYTotal;

    model->accessing_point_data.lock();

    double sim_time = model->prms.SimulationTime;
    const double &h = model->prms.cellsize;
    const double ox = h * model->prms.ModeledRegionOffsetX;
    const double oy = h * model->prms.ModeledRegionOffsetY;

    const int nPts = model->gpu.hssoa.size;
    for(int i=0;i<nPts;i++)
    {
        SOAIterator s = model->gpu.hssoa.begin()+i;
        PointVector2r pos = s->getPos(model->prms.cellsize);
        points->SetPoint((vtkIdType)i, pos[0]+ox, pos[1]+oy, 1.0);
    }
    points->Modified();
    actor_points->GetProperty()->SetPointSize(model->prms.ParticleViewSize);
    double range = std::pow(10,ranges[VisualizingVariable]);
    double centerVal = 0;

    for(int i=0;i<gx;i++)
        for(int j=0;j<gy;j++)
        {
            int idx2 = (j + gy*i);
            if(model->gpu.grid_status_buffer[idx2])
            {
                // water color
                unsigned char* pixel = static_cast<unsigned char*>(uniformGrid->GetScalarPointer(i+_ox, j+_oy, 0));
                pixel[0] = 0x15;
                pixel[1] = 0x1f;
                pixel[2] = 0x2f;
            }
        }


    actor_points->VisibilityOn();

    if(VisualizingVariable == VisOpt::none)
    {
        actor_points->VisibilityOff();
    }
    else if(VisualizingVariable == VisOpt::status)
    {
        scalarBar->VisibilityOn();
        points_mapper->ScalarVisibilityOn();
        points_mapper->SetColorModeToMapScalars();
        points_mapper->UseLookupTableScalarRangeOn();
        for(int i=0;i<nPts;i++)
        {
            float value = 1;
            SOAIterator s = model->gpu.hssoa.begin()+i;
            if(s->getDisabledStatus()) value = 3;
            else if(s->getCrushedStatus()) value = 0;
            else if(s->getWeakenedStatus()) value = 2;
        }
    }

    else if(VisualizingVariable == VisOpt::v_u)
    {
        actor_points->VisibilityOff();
        for(int i=0;i<gx;i++)
            for(int j=0;j<gy;j++)
            {
                int idx2 = (j + gy*i);
                if(model->gpu.grid_status_buffer[idx2])
                {
                    // water color
                    Eigen::Vector2f v = model->fluent_interpolatror.getInterpolation(i,j);
                    float val = (v.x())/range + 0.5;
                    //float val = (sin(i/100.)+1.)/2.;
                    std::array<uint8_t, 3> c = colormap.getColor(ColorMap::Palette::SpecialJ, val);

                    unsigned char* pixel = static_cast<unsigned char*>(uniformGrid->GetScalarPointer(i+_ox, j+_oy, 0));
                    for(int k=0;k<3;k++) pixel[k] = c[k];
                }
            }
    }
    else if(VisualizingVariable == VisOpt::v_v)
    {
        actor_points->VisibilityOff();
        for(int i=0;i<gx;i++)
            for(int j=0;j<gy;j++)
            {
                int idx2 = (j + gy*i);
                if(model->gpu.grid_status_buffer[idx2])
                {
                    // water color
                    Eigen::Vector2f v = model->fluent_interpolatror.getInterpolation(i,j);
                    float val = (v.y())/range + 0.5;
                    //float val = (sin(i/100.)+1.)/2.;
                    std::array<uint8_t, 3> c = colormap.getColor(ColorMap::Palette::SpecialJ, val);

                    unsigned char* pixel = static_cast<unsigned char*>(uniformGrid->GetScalarPointer(i+_ox, j+_oy, 0));
                    for(int k=0;k<3;k++) pixel[k] = c[k];
                }
            }

    }
    else if(VisualizingVariable == VisOpt::v_norm)
    {
        actor_points->VisibilityOff();
        for(int i=0;i<gx;i++)
            for(int j=0;j<gy;j++)
            {
                int idx2 = (j + gy*i);
                if(model->gpu.grid_status_buffer[idx2])
                {
                    // water color
                    Eigen::Vector2f v = model->fluent_interpolatror.getInterpolation(i,j);
                    float val = (v.norm())/range;
                    //float val = (sin(i/100.)+1.)/2.;
                    std::array<uint8_t, 3> c = colormap.getColor(ColorMap::Palette::P2, val);
                    unsigned char* pixel = static_cast<unsigned char*>(uniformGrid->GetScalarPointer(i+_ox, j+_oy, 0));
                    for(int k=0;k<3;k++) pixel[k] = c[k];
                }
            }
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
    }
    else if(VisualizingVariable == VisOpt::Jp_inv)
    {
        scalarBar->VisibilityOn();
        points_mapper->ScalarVisibilityOn();
        points_mapper->SetColorModeToDirectScalars();
        for(int i=0;i<nPts;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            double val = s->getValue(SimParams::idx_Jp_inv)-1.;
            double value = (val)/range + 0.5;
            double alpha = abs(val)/range;

            int pt_idx = s->getValueInt(SimParams::integer_point_idx);
            uint32_t rgb = model->gpu.point_colors_rgb[pt_idx];
            std::array<uint8_t, 3> c = colormap.getColor(ColorMap::Palette::Pressure, value);
            std::array<uint8_t, 3> c2 = colormap.mergeColors(rgb, c, alpha);

            for(int k=0;k<3;k++) pts_colors->SetValue((vtkIdType)(i*3+k), c2[k]);
        }
        points_polydata->GetPointData()->SetActiveScalars("pts_colors");
        pts_colors->Modified();
//        points_polydata->Modified();
        points_filter->Update();
    }
    else if(VisualizingVariable == VisOpt::P)
    {
        scalarBar->VisibilityOn();
        points_mapper->ScalarVisibilityOn();
        points_mapper->SetColorModeToDirectScalars();
        for(int i=0;i<nPts;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            double val = s->getValue(SimParams::idx_P);
            double value = (val)/range+0.5;
            double alpha = abs(val)/range;

            int pt_idx = s->getValueInt(SimParams::integer_point_idx);
            uint32_t rgb = model->gpu.point_colors_rgb[pt_idx];
            std::array<uint8_t, 3> c = colormap.getColor(ColorMap::Palette::Pressure, value);
            std::array<uint8_t, 3> c2 = colormap.mergeColors(rgb, c, alpha);

            for(int k=0;k<3;k++) pts_colors->SetValue((vtkIdType)(i*3+k), c2[k]);
        }
        points_polydata->GetPointData()->SetActiveScalars("pts_colors");
        pts_colors->Modified();
        //        points_polydata->Modified();
        points_filter->Update();
    }
    else if(VisualizingVariable == VisOpt::Q)
    {

    }


    else if(VisualizingVariable == VisOpt::strength)
    {

    }
    else
    {
        points_mapper->ScalarVisibilityOff();
        scalarBar->VisibilityOff();
    }

    points_filter->Update();
    uniformGrid->Modified();
    mapper_uniformgrid->Update();


    model->accessing_point_data.unlock();

}


void icy::VisualRepresentation::ChangeVisualizationOption(int option)
{
    VisualizingVariable = (VisOpt)option;
    SynchronizeTopology();
}




