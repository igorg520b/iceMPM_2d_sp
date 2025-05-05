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
}



void icy::VisualRepresentation::SynchronizeTopology()
{
    LOGV("SynchronizeTopology()");

    const int &width = model->prms.InitializationImageSizeX;
    const int &height = model->prms.InitializationImageSizeY;
    const int &ox = model->prms.ModeledRegionOffsetX;
    const int &oy = model->prms.ModeledRegionOffsetY;
    const int &gx = model->prms.GridXTotal;
    const int &gy = model->prms.GridYTotal;
    const double &h = model->prms.cellsize;

    if(!model->gpu.original_image_colors_rgb.size()) return;

    // (1) background image (but cover the modelled area in blue)
    renderedImage.assign(model->gpu.original_image_colors_rgb.begin(), model->gpu.original_image_colors_rgb.end());

    for(size_t i=0;i<gx;i++)
        for(size_t j=0;j<gy;j++)
        {
            uint8_t status = model->gpu.grid_status_buffer[j + i*gy];
            if(status == 100)
            {
                if(VisualizingVariable == VisOpt::v_norm)
                {
                    // visualize current flow
                    size_t idx2 = j + i*gy;
                    t_GridReal vx = model->wac_interpolator.current_flow_data[idx2];
                    t_GridReal vy = model->wac_interpolator.current_flow_data[idx2 + gx*gy];
                    float norm = sqrt(vx*vx + vy*vy);
                    norm /= 0.3;

                    std::array<uint8_t, 3> c = colormap.getColor(ColorMap::Palette::ANSYS, norm);
                    for(int k=0;k<3;k++) renderedImage[((i+ox) + (j+oy)*width)*3+k] = c[k];
                }
                else if(VisualizingVariable == VisOpt::v_u)
                {
                    // visualize current flow
                    size_t idx2 = j + i*gy;
                    t_GridReal vx = model->wac_interpolator.current_flow_data[idx2];
                    std::array<uint8_t, 3> c = colormap.getColor(ColorMap::Palette::ANSYS, 0.5+vx/0.3);
                    for(int k=0;k<3;k++) renderedImage[((i+ox) + (j+oy)*width)*3+k] = c[k];
                }
                else if(VisualizingVariable == VisOpt::v_v)
                {
                    // visualize current flow
                    size_t idx2 = j + i*gy;
                    t_GridReal vy = model->wac_interpolator.current_flow_data[idx2 + gx*gy];
                    std::array<uint8_t, 3> c = colormap.getColor(ColorMap::Palette::ANSYS, 0.5+vy/0.3);
                    for(int k=0;k<3;k++) renderedImage[((i+ox) + (j+oy)*width)*3+k] = c[k];
                }
                else
                {
//                    for(int k=0;k<3;k++)
//                        renderedImage[((i+ox) + (j+oy)*width)*3+k] = waterColor[k];
                }
            }
            else
            {
                // not modelled area
                if(VisualizingVariable == VisOpt::regions)
                {
                    // visualize current flow
                    size_t idx2 = j + i*gy;
                    uint8_t region_id = model->gpu.grid_status_buffer[idx2];
                    float val = (region_id % 13) / 12.;
                    std::array<uint8_t, 3> c = colormap.getColor(ColorMap::Palette::Pastel, val);
                    for(int k=0;k<3;k++) renderedImage[((i+ox) + (j+oy)*width)*3+k] = c[k];
                }

            }
        }




    raster_scalars->SetNumberOfComponents(3);          // RGB has 3 components
    raster_scalars->SetArray(renderedImage.data(), renderedImage.size(), 1);
    raster_scalars->Modified();

    raster_imageData->SetDimensions(width, height, 1); // 2D image, depth = 1
    raster_imageData->SetSpacing(1.0, 1.0, 1.0);      // Pixel spacing
    raster_imageData->SetOrigin(0.0, 0.0, 0.0);       // Origin at (0,0,0)
    raster_imageData->GetPointData()->SetScalars(raster_scalars);

    raster_plane->SetOrigin(-h/2, -h/2, -1.0);           // Bottom-left corner
    raster_plane->SetPoint1((width-0.5)*h, -h/2, -1.0);         // Bottom-right (x-axis)
    raster_plane->SetPoint2(-h/2, (height-0.5)*h, -1.0);        // Top-left (y-axis)
    raster_plane->SetNormal(0.0, 0.0, 1.0);           // Normal along z-axis (facing forward)

    raster_mapper->SetInputConnection(raster_plane->GetOutputPort());

    raster_texture->SetInputData(raster_imageData);
    raster_texture->InterpolateOff(); // Smooth texture rendering

    raster_actor->SetMapper(raster_mapper);
    raster_actor->SetTexture(raster_texture);

    raster_mapper->Update();
    raster_texture->Update();



    // (2) points
    const int nPts = model->gpu.hssoa.size;
    if(!nPts) return;

    model->accessing_point_data.lock();
    points->SetNumberOfPoints(nPts);
    pts_colors->SetNumberOfValues(nPts*3);

    model->accessing_point_data.unlock();
    SynchronizeValues();
}


void icy::VisualRepresentation::SynchronizeValues()
{
    const int &width = model->prms.InitializationImageSizeX;
    const int &height = model->prms.InitializationImageSizeY;
    const int &ox = model->prms.ModeledRegionOffsetX;
    const int &oy = model->prms.ModeledRegionOffsetY;
    const int &gx = model->prms.GridXTotal;
    const int &gy = model->prms.GridYTotal;
    const double &h = model->prms.cellsize;

    model->accessing_point_data.lock();

    // points' coordinates
    const int nPts = model->gpu.hssoa.size;
    for(int i=0;i<nPts;i++)
    {
        SOAIterator s = model->gpu.hssoa.begin()+i;
        PointVector2r pos = s->getPos(model->prms.cellsize);
        points->SetPoint((vtkIdType)i, pos[0]+ox*h, pos[1]+oy*h, 1.0);
    }


    // point visibility and color
    actor_points->VisibilityOn();
    scalarBar->VisibilityOff();
    points_mapper->ScalarVisibilityOn();
    points_mapper->SetColorModeToDirectScalars();
    points_mapper->Modified();

    const double range = std::pow(10,ranges[VisualizingVariable]);

    if(VisualizingVariable == VisOpt::color)
    {
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
    }
    else if(VisualizingVariable == VisOpt::status)
    {
        for(int i=0;i<nPts;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            bool crushed = s->getCrushedStatus();
            bool disabled = s->getDisabledStatus();

            double val = 0;
            if(crushed) val = 1;
            if(disabled) val = 2;
            std::array<uint8_t, 3> c = colormap.getColor(ColorMap::Palette::NCD, val/2.);
            for(int k=0;k<3;k++) pts_colors->SetValue((vtkIdType)(i*3+k), c[k]);
        }
    }
    else if(VisualizingVariable == VisOpt::none || VisualizingVariable == VisOpt::regions)
    {
        for(int i=0;i<nPts;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            pts_colors->SetValue((vtkIdType)(i*3+0), 240);
            pts_colors->SetValue((vtkIdType)(i*3+1), 122);
            pts_colors->SetValue((vtkIdType)(i*3+2), 122);
        }
    }
    else if(VisualizingVariable == VisOpt::Jp_inv)
    {
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
    }
    else if(VisualizingVariable == VisOpt::P)
    {
        for(int i=0;i<nPts;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;

            const double val = s->getValue(SimParams::idx_P);
            double value = (val)/range+0.5;
            double alpha = abs(val)/range;

            int pt_idx = s->getValueInt(SimParams::integer_point_idx);
            uint32_t rgb = model->gpu.point_colors_rgb[pt_idx];
            std::array<uint8_t, 3> c = colormap.getColor(ColorMap::Palette::Pressure, value);
            std::array<uint8_t, 3> c2 = colormap.mergeColors(rgb, c, alpha);
            for(int k=0;k<3;k++) pts_colors->SetValue((vtkIdType)(i*3+k), c2[k]);
        }
    }
    else if(VisualizingVariable == VisOpt::Q)
    {
        for(int i=0;i<nPts;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            double val = s->getValue(SimParams::idx_Q);
            double value = (val)/range+0.5;
            double alpha = abs(val)/range;

            int pt_idx = s->getValueInt(SimParams::integer_point_idx);
            uint32_t rgb = model->gpu.point_colors_rgb[pt_idx];
            std::array<uint8_t, 3> c = colormap.getColor(ColorMap::Palette::P2, value);
            std::array<uint8_t, 3> c2 = colormap.mergeColors(rgb, c, alpha);

            for(int k=0;k<3;k++) pts_colors->SetValue((vtkIdType)(i*3+k), c2[k]);
        }
    }
    else if(VisualizingVariable == VisOpt::thickness)
    {
        for(int i=0;i<nPts;i++)
        {
            SOAIterator s = model->gpu.hssoa.begin()+i;
            double val = s->getValue(SimParams::idx_thickness);
            double value = (val-0.8)/0.2;
            std::array<uint8_t, 3> c = colormap.getColor(ColorMap::Palette::ANSYS, value);
            for(int k=0;k<3;k++) pts_colors->SetValue((vtkIdType)(i*3+k), c[k]);
        }
    }
    else
    {
        actor_points->VisibilityOff();
    }

    points_filter->Update();

    points_polydata->GetPointData()->SetActiveScalars("pts_colors");
    pts_colors->Modified();
    points_polydata->Modified();
    points->Modified();

    actor_points->GetProperty()->SetPointSize(model->prms.ParticleViewSize);
    model->accessing_point_data.unlock();
}


void icy::VisualRepresentation::ChangeVisualizationOption(int option)
{
    VisualizingVariable = (VisOpt)option;
    SynchronizeTopology();
}




