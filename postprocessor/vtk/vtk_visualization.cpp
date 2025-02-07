#include "vtk_visualization.h"
#include "parameters_sim.h"
//#include <omp.h>
#include <algorithm>
#include <iostream>
#include <spdlog/spdlog.h>
#include <Eigen/Core>


VTKVisualization::VTKVisualization(FrameData& fd) : frameData(fd)
{
    colormap.populateLut(ColorMap::Palette::Pressure, lut_Pressure);

    // text
    constexpr int fontSize = 20;
    vtkTextProperty* txtprop = actor_text->GetTextProperty();
    txtprop->SetFontFamilyToArial();
    txtprop->BoldOff();
    txtprop->SetFontSize(fontSize);
    txtprop->ShadowOff();
    txtprop->SetColor(0,0,0);

    // main grid
    mapper_uniformgrid->SetInputData(uniformGrid);
    actor_grid_main->SetMapper(mapper_uniformgrid);
    actor_grid_main->GetProperty()->SetEdgeVisibility(false);
    actor_grid_main->GetProperty()->LightingOff();
    actor_grid_main->GetProperty()->ShadingOff();
    actor_grid_main->PickableOff();


    // scalar bar
    scalarBar->SetLookupTable(lut_Pressure);
    scalarBar->SetMaximumWidthInPixels(180);
    scalarBar->SetBarRatio(0.1);
    scalarBar->SetMaximumHeightInPixels(250);
    scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();
    scalarBar->GetPositionCoordinate()->SetValue(0.01,0.015, 0.0);
    scalarBar->SetLabelFormat("%.1e");
    scalarBar->GetLabelTextProperty()->BoldOff();
    scalarBar->GetLabelTextProperty()->ItalicOff();
    scalarBar->GetLabelTextProperty()->ShadowOff();
    scalarBar->GetLabelTextProperty()->SetColor(0.1,0.1,0.1);
    scalarBar->GetLabelTextProperty()->SetFontFamilyToTimes();
    scalarBar->SetUnconstrainedFontSize(true);
    scalarBar->GetLabelTextProperty()->SetFontSize(40);


    // rectangle frame
    rectangleActor->SetMapper(rectangleMapper);
    rectangleActor->GetProperty()->SetColor(0.0, 0.0, 0.0); // Red color
    rectangleActor->GetProperty()->SetLineWidth(2.0);       // Line thickness

    rectangleMapper->SetInputData(rectanglePolyData);
    rectanglePolyData->SetPoints(rectanglePoints);
    rectanglePolyData->SetLines(rectangleLines);

    actor_text->SetDisplayPosition(600, 10);

    txtprop = actor_text_title->GetTextProperty();
    txtprop->SetFontFamilyToArial();
    txtprop->BoldOff();
    txtprop->ShadowOff();
    txtprop->SetColor(0,0,0);
    txtprop->SetFontSize(fontSize);

    actor_text_title->SetDisplayPosition(10, 500);

}



void VTKVisualization::SynchronizeTopology()
{
    spdlog::info("VTKVisualization::SynchronizeTopology()");
    if(!frameData.dataLoaded) return;

    // actor_grid_main
    const int gx = frameData.prms.GridXTotal;
    const int gy = frameData.prms.GridYTotal;
    const double h = frameData.prms.cellsize;

    const int imgx = frameData.prms.InitializationImageSizeX;
    const int imgy = frameData.prms.InitializationImageSizeY;

    int ox = frameData.prms.ModeledRegionOffsetX;
    int oy = frameData.prms.ModeledRegionOffsetY;

    spdlog::info("img size {}", frameData.original_image_colors_rgb.size());
    spdlog::info("calculated size {}", imgx*imgy*3);

    uniformGrid->SetDimensions(imgx, imgy, 1);
    uniformGrid->SetSpacing(h, h, 1.0);
    uniformGrid->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    for(int i=0;i<imgx;i++)
        for(int j=0;j<imgy;j++)
        {
            unsigned char* pixel = static_cast<unsigned char*>(uniformGrid->GetScalarPointer(i, j, 0));
            int idx = (i+imgx*j)*3;
            pixel[0] = frameData.original_image_colors_rgb[idx+0];
            pixel[1] = frameData.original_image_colors_rgb[idx+1];
            pixel[2] = frameData.original_image_colors_rgb[idx+2];
        }

    for(int i=0;i<gx;i++)
        for(int j=0;j<gy;j++)
        {
            int idx2 = (j + gy*i);
            if(frameData.grid_status_buffer[idx2])
            {
                // water color
                unsigned char* pixel = static_cast<unsigned char*>(uniformGrid->GetScalarPointer(i+ox, j+oy, 0));
                pixel[0] = 0x15;
                pixel[1] = 0x1f;
                pixel[2] = 0x2f;
            }
        }

    mapper_uniformgrid->Update();


    // rectangle
    rectanglePoints->SetNumberOfPoints(4);
    rectanglePoints->SetPoint(0, 0.5*h, 0.5*h,2.0);
    rectanglePoints->SetPoint(1, (imgx-1.5)*h, 0.5*h , 2.0);
    rectanglePoints->SetPoint(2, (imgx-1.5)*h, (imgy-1.5)*h , 2.0);
    rectanglePoints->SetPoint(3, 0.5*h, (imgy-1.5)*h , 2.0);

    vtkIdType pointIds[5] = {0, 1, 2, 3, 0};
    rectangleLines->Reset();
    rectangleLines->InsertNextCell(5, pointIds);


    SynchronizeValues();
}


void VTKVisualization::SynchronizeValues()
{
    if(!frameData.dataLoaded) return;
    const int &_ox = frameData.prms.ModeledRegionOffsetX;
    const int &_oy = frameData.prms.ModeledRegionOffsetY;

    const int &gx = frameData.prms.GridXTotal;
    const int &gy = frameData.prms.GridYTotal;

    double sim_time = frameData.prms.SimulationTime;
    const double &h = frameData.prms.cellsize;
    const double ox = h * frameData.prms.ModeledRegionOffsetX;
    const double oy = h * frameData.prms.ModeledRegionOffsetY;

    double range = std::pow(10,ranges[VisualizingVariable]);
    double centerVal = 0;
    spdlog::info("SynchronizeValues() {}; range {}", this->VisualizingVariable, range);

    const std::array<uint8_t, 3>& seaWater {0x15, 0x1f, 0x2f};

    for(int i=0;i<gx;i++)
        for(int j=0;j<gy;j++)
        {
            int idx2 = (j + gy*i);
            if(frameData.grid_status_buffer[idx2])
            {
                // water color
                unsigned char* pixel = static_cast<unsigned char*>(uniformGrid->GetScalarPointer(i+_ox, j+_oy, 0));
                pixel[0] = 0x15;
                pixel[1] = 0x1f;
                pixel[2] = 0x2f;
            }
        }


    if(VisualizingVariable == VisOpt::Jp_inv)
    {
        spdlog::info("rendering Jp_inv");
        for(int i=0;i<gx;i++)
            for(int j=0;j<gy;j++)
            {
                int idx2 = (j + gy*i);
                if(frameData.grid_status_buffer[idx2])
                {
                    uint8_t count = frameData.count[idx2];
                    float Jp_inv = frameData.vis_Jpinv[idx2];
                    std::array<uint8_t, 3> _rgb;
                    for(int k=0;k<3;k++) _rgb[k] = frameData.rgb[idx2*3+k];

                    if(count > 0)
                    {
                        float coeff2 = std::clamp(std::abs(Jp_inv-1)/range, 0., 1.);
                        coeff2 = pow(coeff2, 2);
                        float coeff1 = std::min(count/2.,1.); // how much water surface is obscured
                        float val = (Jp_inv-1.)/range + 0.5;
                        if(Jp_inv>1.) coeff2 = std::clamp(std::abs(Jp_inv-1)*0.5/range, 0., 1.);
                        std::array<uint8_t, 3> c = colormap.getColor(ColorMap::Palette::Pressure, val);

                        std::array<uint8_t, 3> c2 = ColorMap::mergeColors(_rgb, c, coeff2);
                        std::array<uint8_t, 3> c3 = ColorMap::mergeColors(seaWater, c2, coeff1);

                        unsigned char* pixel = static_cast<unsigned char*>(uniformGrid->GetScalarPointer(i+_ox, j+_oy, 0));
                        for(int k=0;k<3;k++) pixel[k] = c3[k];
                    }
                }
            }

        scalarBar->VisibilityOn();
        scalarBar->SetLookupTable(lut_Pressure);
        scalarBar->SetLabelFormat("%.1f");
        lut_Pressure->SetTableRange(1-range, 1+range);

    }
    else if(VisualizingVariable == VisOpt::colors)
    {
        spdlog::info("rendering colors");

        for(int i=0;i<gx;i++)
            for(int j=0;j<gy;j++)
            {
                int idx2 = (j + gy*i);
                if(frameData.grid_status_buffer[idx2])
                {
                    uint8_t count = frameData.count[idx2];
                    std::array<uint8_t, 3> _rgb;
                    for(int k=0;k<3;k++) _rgb[k] = frameData.rgb[idx2*3+k];

                    if(count > 0)
                    {
                        unsigned char* pixel = static_cast<unsigned char*>(uniformGrid->GetScalarPointer(i+_ox, j+_oy, 0));
                        for(int k=0;k<3;k++) pixel[k] = _rgb[k];
                    }
                }
            }
    }
    uniformGrid->Modified();
        mapper_uniformgrid->Update();

    int64_t display_date = (int64_t)frameData.prms.SimulationTime + frameData.prms.SimulationStartUnixTime;
    std::time_t unix_time = display_date;
    std::tm* tm_time = std::gmtime(&unix_time);
    // Format the time
    char buffer[100];
    std::strftime(buffer, sizeof(buffer), "%d-%m-%Y %H:%M:%S UTC", tm_time);
    actor_text->SetInput(buffer);


    /*
    std::vector<uint8_t> &gb = frameData.grid_status_buffer;
    if(!frameData.dataLoaded || gb.size()==0) return;

    const int gx = frameData.prms.GridXTotal;
    const int gy = frameData.prms.GridY;
    const double range = std::pow(10,ranges[VisualizingVariable]); // visualization scale

    const Eigen::Vector3f land_color(0.35, 0.15, 0.15);
    const Eigen::Vector3f water_color(0x11/255., 0x18/255., 0x20/255.);


    const double h = frameData.prms.cellsize;
    const float spacingX = gx*h*0.9/(numwCols-1);
    const float spacingY = gy*h*0.9/(numwRows-1);
    const float wind_scale = gx*h*0.002;

    frameData.windInterpolator.prepareVisualizationData(frameData.prms.SimulationTime);
    for (int i = 0; i < numwCols; ++i)
    {
        for (int j = 0; j < numwRows; ++j)
        {
            int idx = i+j*numwCols;
            float x = i * spacingX + gx*h*0.05;
            float y = j * spacingY + gy*h*0.05;

            float lat = frameData.prms.LatMin + (frameData.prms.LatMax-frameData.prms.LatMin)*y/(gy*h);
            float lon = frameData.prms.LonMin + (frameData.prms.LonMax-frameData.prms.LonMin)*x/(gx*h);

            Eigen::Vector2f wv = frameData.windInterpolator.vis_interp(lat, lon);
            //wv*=wind_scale;

            vectors_values->SetValue(idx*3+0, wv.x());
            vectors_values->SetValue(idx*3+1, wv.y());
        }
    }
    vectors_values->Modified();




    if(VisualizingVariable == VisOpt::none)
    {
        scalarBar->VisibilityOff();
        for(int idx_y=0; idx_y<gy; idx_y++)
            for(int idx_x=0; idx_x<gx; idx_x++)
            {
                int idx_visual = idx_x + idx_y*gx;
                int idx_storage = idx_y + idx_x*gy;
                int value = gb[idx_storage];
                visualized_values_grid->SetValue(idx_visual, value);
            }
        visualized_values_grid->Modified();

        mapper_grid_main->SetLookupTable(hueLut_count);
        mapper_grid_main->SetColorModeToMapScalars();
        visualized_values_grid->Modified();
        structuredGrid_main->GetCellData()->SetActiveScalars("visualized_values_grid");

    }
    else if(VisualizingVariable == VisOpt::count)
    {
        scalarBar->VisibilityOff();
        for(int idx_y=0; idx_y<gy; idx_y++)
            for(int idx_x=0; idx_x<gx; idx_x++)
            {
                int idx_visual = idx_x + idx_y*gx;
                int idx_storage = idx_y + idx_x*gy;
                int value = gb[idx_storage];
                if(value == 0 && frameData.count[idx_storage]>0) value = 1 + frameData.count[idx_storage];
                visualized_values_grid->SetValue(idx_visual, value);
            }
        visualized_values_grid->Modified();

        mapper_grid_main->SetLookupTable(hueLut_count);
        mapper_grid_main->SetColorModeToMapScalars();
        structuredGrid_main->GetCellData()->SetActiveScalars("visualized_values_grid");

        visualized_values_grid->Modified();
    }
    else if(VisualizingVariable == VisOpt::colors)
    {
        scalarBar->VisibilityOff();
        for(int idx_y=0; idx_y<gy; idx_y++)
            for(int idx_x=0; idx_x<gx; idx_x++)
            {
                int idx_visual = idx_x + idx_y*gx;
                int idx_storage = idx_y + idx_x*gy;
                int land = gb[idx_storage];
                int count = frameData.count[idx_storage];
                float r = frameData.vis_r[idx_storage];
                float g = frameData.vis_g[idx_storage];
                float b = frameData.vis_b[idx_storage];

                float coeff = 0;
                if(!land) coeff = std::min(count/3.,1.);

                if(land)
                {
                    pts_colors->SetValue((vtkIdType)(idx_visual*3+0), 0.35*255);
                    pts_colors->SetValue((vtkIdType)(idx_visual*3+1), 0.15*255);
                    pts_colors->SetValue((vtkIdType)(idx_visual*3+2), 0.15*255);
                }
                else
                {
                    unsigned char ur = (unsigned char)((1-coeff)*0x11 + coeff*255*r);
                    unsigned char ug = (unsigned char)((1-coeff)*0x18 + coeff*255*g);
                    unsigned char ub = (unsigned char)((1-coeff)*0x20 + coeff*255*b);
                    pts_colors->SetValue((vtkIdType)(idx_visual*3+0), ur);
                    pts_colors->SetValue((vtkIdType)(idx_visual*3+1), ug);
                    pts_colors->SetValue((vtkIdType)(idx_visual*3+2), ub);
                }
            }
        visualized_values_grid->Modified();
        structuredGrid_main->GetCellData()->SetActiveScalars("pts_colors");
        mapper_grid_main->SetColorModeToDirectScalars();
    }
    else if(VisualizingVariable == VisOpt::Jp_inv)
    {
        for(int idx_y=0; idx_y<gy; idx_y++)
            for(int idx_x=0; idx_x<gx; idx_x++)
            {
                int idx_visual = idx_x + idx_y*gx;
                int idx_storage = idx_y + idx_x*gy;
                int land = gb[idx_storage];
                int count = frameData.count[idx_storage];
                float r = frameData.vis_r[idx_storage];
                float g = frameData.vis_g[idx_storage];
                float b = frameData.vis_b[idx_storage];
                //float a = frameData.vis_alpha[idx_storage];
                float Jp_inv = frameData.vis_Jpinv[idx_storage];


                const Eigen::Vector3f land_color(0.35, 0.15, 0.15);
                const Eigen::Vector3f water_color(0x11/255., 0x18/255., 0x20/255.);
                const Eigen::Vector3f ice_cover_color(r,g,b);
                Eigen::Vector3f paint_color(0,0,0);



                if(land)
                {
                    paint_color = land_color;
                }
                else
                {
                    float coeff1 = std::min(count/3.,1.); // how much water surface is obscured
                    paint_color = (1-coeff1)*water_color + coeff1*ice_cover_color;

                    if(!std::isnan(Jp_inv))
                    {
                        float coeff2 = 0.05+std::abs(Jp_inv-1)/range;
                        coeff2 = std::clamp(coeff2, 0.25f, 1.0f);
                        Eigen::Vector3f Jpinv_color = interpolateColor(lutSpecialJ, 5, 0.5+(Jp_inv-1)/range);
                        paint_color = (1-coeff2)*paint_color + coeff2*Jpinv_color;
                    }


                }

                for(int k=0;k<3;k++) pts_colors->SetValue((vtkIdType)(idx_visual*3+k), (unsigned char)(255*paint_color[k]));
            }
        visualized_values_grid->Modified();
        structuredGrid_main->GetCellData()->SetActiveScalars("pts_colors");
        mapper_grid_main->SetColorModeToDirectScalars();

        scalarBar->VisibilityOn();
        scalarBar->SetLookupTable(hueLut_J);
        scalarBar->SetLabelFormat("%.1f");
        hueLut_J->SetTableRange(1-range, 1+range);
    }
    else if(VisualizingVariable == VisOpt::P)
    {
        for(int idx_y=0; idx_y<gy; idx_y++)
            for(int idx_x=0; idx_x<gx; idx_x++)
            {
                int idx_visual = idx_x + idx_y*gx;
                int idx_storage = idx_y + idx_x*gy;
                int land = gb[idx_storage];
                int count = frameData.count[idx_storage];
                float r = frameData.vis_r[idx_storage];
                float g = frameData.vis_g[idx_storage];
                float b = frameData.vis_b[idx_storage];
                float P = frameData.vis_P[idx_storage];

                const Eigen::Vector3f ice_cover_color(r,g,b);
                Eigen::Vector3f paint_color(0,0,0);

                if(land)
                {
                    paint_color = land_color;
                }
                else
                {
                    float coeff1 = std::min((double)count/3.,1.); // how much water surface is obscured
                    paint_color = (1-coeff1)*water_color + coeff1*ice_cover_color;

                    if(!std::isnan(P))
                    {
                        //float coeff2 = 0.00+std::abs(P)/range;
                        float coeff2 = 0.00+P*P/(range*range*3);
                        coeff2 = std::clamp(coeff2, 0.f, 1.0f);
                        Eigen::Vector3f P_color = interpolateColor(naturalRidges, 7, 0.5+(P)/range);
                        paint_color = (1-coeff2)*paint_color + coeff2*P_color;
                    }
                }
                for(int k=0;k<3;k++) pts_colors->SetValue((vtkIdType)(idx_visual*3+k), (unsigned char)(255*paint_color[k]));
            }
        visualized_values_grid->Modified();
        structuredGrid_main->GetCellData()->SetActiveScalars("pts_colors");
        mapper_grid_main->SetColorModeToDirectScalars();

        scalarBar->VisibilityOn();
        scalarBar->SetLookupTable(hueLut_Jpinv);
        scalarBar->SetLabelFormat("%.1e");
        hueLut_Jpinv->SetTableRange(-range, range);
    }

    else if(VisualizingVariable == VisOpt::Q)
    {
        for(int idx_y=0; idx_y<gy; idx_y++)
            for(int idx_x=0; idx_x<gx; idx_x++)
            {
                int idx_visual = idx_x + idx_y*gx;
                int idx_storage = idx_y + idx_x*gy;
                int land = gb[idx_storage];
                int count = frameData.count[idx_storage];
                float r = frameData.vis_r[idx_storage];
                float g = frameData.vis_g[idx_storage];
                float b = frameData.vis_b[idx_storage];
                float Q = frameData.vis_Q[idx_storage];

                const Eigen::Vector3f ice_cover_color(r,g,b);
                Eigen::Vector3f paint_color(0,0,0);

                if(land)
                {
                    paint_color = land_color;
                }
                else
                {
                    float coeff1 = std::min((double)count/3.,1.); // how much water surface is obscured
                    paint_color = (1-coeff1)*water_color + coeff1*ice_cover_color;

                    if(!std::isnan(Q))
                    {
                        //float coeff2 = 0.00+std::abs(P)/range;
                        float coeff2 = 0.00+Q*Q/(range*range*3);
                        coeff2 = std::clamp(coeff2, 0.f, 1.0f);
                        Eigen::Vector3f Q_color = interpolateColor(lutSpecialQ, 9, Q/range);
                        paint_color = (1-coeff2)*paint_color + coeff2*Q_color;
                    }
                }
                for(int k=0;k<3;k++) pts_colors->SetValue((vtkIdType)(idx_visual*3+k), (unsigned char)(255*paint_color[k]));
            }
        visualized_values_grid->Modified();
        structuredGrid_main->GetCellData()->SetActiveScalars("pts_colors");
        mapper_grid_main->SetColorModeToDirectScalars();

        scalarBar->VisibilityOn();
        scalarBar->SetLookupTable(hueLut_Q);
        scalarBar->SetLabelFormat("%.1e");
        hueLut_Q->SetTableRange(0, range);
    }
    else if(VisualizingVariable == VisOpt::wind_norm)
    {
        //frameData.windInterpolator.setTime(frameData.prms.SimulationTime);
        //float tb = frameData.windInterpolator.interpolationCoeffFromTime(frameData.prms.SimulationTime);

        frameData.windInterpolator.prepareVisualizationData(frameData.prms.SimulationTime);

        for(int idx_y=0; idx_y<gy; idx_y++)
            for(int idx_x=0; idx_x<gx; idx_x++)
            {
                int idx_visual = idx_x + idx_y*gx;

                float lat = frameData.prms.LatMin + (frameData.prms.LatMax-frameData.prms.LatMin)*(float)idx_y/(float)gy;
                float lon = frameData.prms.LonMin + (frameData.prms.LonMax-frameData.prms.LonMin)*(float)idx_x/(float)gx;

                //Eigen::Vector2f wv = frameData.windInterpolator.interpolationResult(lat, lon, tb);

                Eigen::Vector2f wv = frameData.windInterpolator.vis_interp(lat, lon);
                visualized_values_grid->SetValue(idx_visual, wv.norm());
            }
        mapper_grid_main->SetLookupTable(hueLut);
        mapper_grid_main->SetColorModeToMapScalars();
        structuredGrid_main->GetCellData()->SetActiveScalars("visualized_values_grid");

        scalarBar->VisibilityOn();
        scalarBar->SetLookupTable(hueLut);
        scalarBar->SetLabelFormat("%.1f");
        hueLut->SetTableRange(0, range);
    }





*/
}


void VTKVisualization::ChangeVisualizationOption(int option)
{
    VisualizingVariable = (VisOpt)option;
    if(option == VisOpt::P)
    {
        actor_text_title->SetInput("Hydrostatic pressure");
    }
    else if(option == VisOpt::Q)
    {
        actor_text_title->SetInput("Deviatoric stress");
    }
    else if(option == VisOpt::Jp_inv)
    {
        actor_text_title->SetInput("Rel.surf.dens.");
    }
    else if(option == VisOpt::wind_norm)
    {
        actor_text_title->SetInput("Wind speed");
    }
    else
    {
        actor_text_title->SetInput("-");
    }

    SynchronizeValues();
}


