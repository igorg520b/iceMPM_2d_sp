#include "vtk_visualization.h"
#include "parameters_sim.h"
//#include <omp.h>
#include <algorithm>
#include <iostream>
#include <spdlog/spdlog.h>
#include <Eigen/Core>


VTKVisualization::VTKVisualization(FrameData& fd) : frameData(fd)
{
    populateLut(lutSpecialCount, (sizeof(lutSpecialCount)/(sizeof(lutSpecialCount[0]))), hueLut_count);
    interpolateLut(lutSpecialJ, (sizeof(lutSpecialJ)/(sizeof(lutSpecialJ[0]))), hueLut_J);
    interpolateLut(naturalRidges, (sizeof(naturalRidges)/(sizeof(naturalRidges[0]))), hueLut_Jpinv);
    interpolateLut(lutSpecialQ, (sizeof(lutSpecialQ)/(sizeof(lutSpecialQ[0]))), hueLut_Q);

    // text
    constexpr int fontSize = 20;
    vtkTextProperty* txtprop = actor_text->GetTextProperty();
    txtprop->SetFontFamilyToArial();
    txtprop->BoldOff();
    txtprop->SetFontSize(fontSize);
    txtprop->ShadowOff();
    txtprop->SetColor(0,0,0);

    // main grid
    mapper_grid_main->SetInputData(structuredGrid_main);
    actor_grid_main->SetMapper(mapper_grid_main);

    actor_grid_main->GetProperty()->SetEdgeVisibility(false);
    actor_grid_main->GetProperty()->LightingOff();
    actor_grid_main->GetProperty()->ShadingOff();
    actor_grid_main->PickableOff();

    // main grid color data
    visualized_values_grid->SetName("visualized_values_grid");
    structuredGrid_main->GetCellData()->AddArray(visualized_values_grid);
    mapper_grid_main->UseLookupTableScalarRangeOn();

    pts_colors->SetName("pts_colors");
    pts_colors->SetNumberOfComponents(3);
    structuredGrid_main->GetCellData()->AddArray(pts_colors);

    structuredGrid_main->GetCellData()->SetActiveScalars("visualized_values_grid");


    // scalar bar
    scalarBar->SetLookupTable(hueLut_count);
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

    rectangleMapper_copy1->SetInputData(rectanglePolyData);

    // copies for offscreen renderer
    mapper_grid_main_copy1->SetInputData(structuredGrid_main);
    actor_grid_main_copy1->SetMapper(mapper_grid_main_copy1);
    actor_grid_main_copy1->GetProperty()->SetEdgeVisibility(false);
    actor_grid_main_copy1->GetProperty()->LightingOff();
    actor_grid_main_copy1->GetProperty()->ShadingOff();
    actor_grid_main_copy1->PickableOff();

    vtkTextProperty* txtprop1 = actor_text_copy1->GetTextProperty();
    txtprop1->SetFontFamilyToArial();
    txtprop1->BoldOff();
    txtprop1->ShadowOff();
    txtprop1->SetColor(0,0,0);
    txtprop1->SetFontSize(fontSize);
    actor_text_copy1->SetDisplayPosition(600, 10);
    actor_text->SetDisplayPosition(600, 10);

    txtprop1 = actor_text_title->GetTextProperty();
    txtprop1->SetFontFamilyToArial();
    txtprop1->BoldOff();
    txtprop1->ShadowOff();
    txtprop1->SetColor(0,0,0);
    txtprop1->SetFontSize(fontSize);

    txtprop1 = actor_text_title_copy1->GetTextProperty();
    txtprop1->SetFontFamilyToArial();
    txtprop1->BoldOff();
    txtprop1->ShadowOff();
    txtprop1->SetColor(0,0,0);
    txtprop1->SetFontSize(fontSize);

    actor_text_title->SetDisplayPosition(10, 500);
    actor_text_title_copy1->SetDisplayPosition(10, 500);



    scalarBar_copy1->SetLookupTable(hueLut_count);
    scalarBar_copy1->SetMaximumWidthInPixels(180);
    scalarBar_copy1->SetBarRatio(0.1);
    scalarBar_copy1->SetMaximumHeightInPixels(250);
    scalarBar_copy1->GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();
    scalarBar_copy1->GetPositionCoordinate()->SetValue(0.01,0.015, 0.0);
    scalarBar_copy1->SetLabelFormat("%.1e");
    scalarBar_copy1->GetLabelTextProperty()->BoldOff();
    scalarBar_copy1->GetLabelTextProperty()->ItalicOff();
    scalarBar_copy1->GetLabelTextProperty()->ShadowOff();
    scalarBar_copy1->GetLabelTextProperty()->SetColor(0.1,0.1,0.1);
    scalarBar_copy1->GetLabelTextProperty()->SetFontFamilyToTimes();
    scalarBar_copy1->SetUnconstrainedFontSize(true);
    scalarBar_copy1->GetLabelTextProperty()->SetFontSize(40);

    rectangleActor_copy1->SetMapper(rectangleMapper_copy1);
    rectangleActor_copy1->GetProperty()->SetColor(0.0, 0.0, 0.0); // Red color
    rectangleActor_copy1->GetProperty()->SetLineWidth(2.0);       // Line thickness
}



void VTKVisualization::SynchronizeTopology()
{
    spdlog::info("SynchronizeTopology()");
    if(!frameData.dataLoaded) return;

    // actor_grid_main
    const int gx = frameData.prms.GridXTotal;
    const int gy = frameData.prms.GridY;
    const double h = frameData.prms.cellsize;

    int gx1 = gx+1;
    int gy1 = gy+1;

    structuredGrid_main->SetDimensions(gx1, gy1, 1);
    grid_main_points->SetNumberOfPoints(gx1*gy1);
    for(int idx_y=0; idx_y<=gy; idx_y++)
        for(int idx_x=0; idx_x<=gx; idx_x++)
        {
            double x = (idx_x-0.5)*h;
            double y = (idx_y-0.5)*h;
            double pt_pos[3] {x, y, -3.0};
            grid_main_points->SetPoint((vtkIdType)(idx_x+idx_y*gx1), pt_pos);
        }
    structuredGrid_main->SetPoints(grid_main_points);

    visualized_values_grid->SetNumberOfValues(gx*gy);
    pts_colors->SetNumberOfValues(gx*gy*3);


    // rectangle
    rectanglePoints->SetNumberOfPoints(4);
    rectanglePoints->SetPoint(0, 0.5*h, 0.5*h,-2.0);
    rectanglePoints->SetPoint(1, (gx-1.5)*h, 0.5*h , -2.0);
    rectanglePoints->SetPoint(2, (gx-1.5)*h, (gy-1.5)*h , -2.0);
    rectanglePoints->SetPoint(3, 0.5*h, (gy-1.5)*h , -2.0);

    vtkIdType pointIds[5] = {0, 1, 2, 3, 0};
    rectangleLines->Reset();
    rectangleLines->InsertNextCell(5, pointIds);

    SynchronizeValues();
}


void VTKVisualization::SynchronizeValues()
{
    std::vector<uint8_t> &gb = frameData.grid_status_buffer;
    if(!frameData.dataLoaded || gb.size()==0) return;

    const int gx = frameData.prms.GridXTotal;
    const int gy = frameData.prms.GridY;
    const double range = std::pow(10,ranges[VisualizingVariable]); // visualization scale

    const Eigen::Vector3f land_color(0.35, 0.15, 0.15);
    const Eigen::Vector3f water_color(0x11/255., 0x18/255., 0x20/255.);

    if(VisualizingVariable == VisOpt::none)
    {
        scalarBar->VisibilityOff();
        scalarBar_copy1->VisibilityOff();
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
        scalarBar_copy1->VisibilityOff();
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
        scalarBar_copy1->VisibilityOff();
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
        scalarBar_copy1->VisibilityOn();
        scalarBar_copy1->SetLookupTable(hueLut_J);
        scalarBar_copy1->SetLabelFormat("%.1f");
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
        scalarBar_copy1->VisibilityOn();
        scalarBar_copy1->SetLookupTable(hueLut_Jpinv);
        scalarBar_copy1->SetLabelFormat("%.1e");
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
        scalarBar_copy1->VisibilityOn();
        scalarBar_copy1->SetLookupTable(hueLut_Q);
        scalarBar_copy1->SetLabelFormat("%.1e");
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
        scalarBar_copy1->VisibilityOn();
        scalarBar_copy1->SetLookupTable(hueLut);
        scalarBar_copy1->SetLabelFormat("%.1f");
        hueLut->SetTableRange(0, range);
    }
    else if(VisualizingVariable == VisOpt::divergence)
    {
        frameData.windInterpolator.prepareVisualizationData(frameData.prms.SimulationTime);
        frameData.windInterpolator.computePotentialApproximation();

        for(int idx_y=0; idx_y<gy; idx_y++)
            for(int idx_x=0; idx_x<gx; idx_x++)
            {
                int idx_visual = idx_x + idx_y*gx;

                float lat = frameData.prms.LatMin + (frameData.prms.LatMax-frameData.prms.LatMin)*(float)idx_y/(float)gy;
                float lon = frameData.prms.LonMin + (frameData.prms.LonMax-frameData.prms.LonMin)*(float)idx_x/(float)gx;

                float val = frameData.windInterpolator.vis_interp_divergence(lat, lon);
                visualized_values_grid->SetValue(idx_visual, val);
            }
        mapper_grid_main->SetLookupTable(hueLut);
        mapper_grid_main->SetColorModeToMapScalars();
        structuredGrid_main->GetCellData()->SetActiveScalars("visualized_values_grid");

        scalarBar->VisibilityOn();
        scalarBar->SetLookupTable(hueLut);
        scalarBar->SetLabelFormat("%.1f");
        scalarBar_copy1->VisibilityOn();
        scalarBar_copy1->SetLookupTable(hueLut);
        scalarBar_copy1->SetLabelFormat("%.1f");
        hueLut->SetTableRange(0, range);
    }
    else if(VisualizingVariable == VisOpt::phi)
    {
        frameData.windInterpolator.prepareVisualizationData(frameData.prms.SimulationTime);
        frameData.windInterpolator.computePotentialApproximation();

        for(int idx_y=0; idx_y<gy; idx_y++)
            for(int idx_x=0; idx_x<gx; idx_x++)
            {
                int idx_visual = idx_x + idx_y*gx;

                float lat = frameData.prms.LatMin + (frameData.prms.LatMax-frameData.prms.LatMin)*(float)idx_y/(float)gy;
                float lon = frameData.prms.LonMin + (frameData.prms.LonMax-frameData.prms.LonMin)*(float)idx_x/(float)gx;

                float val = frameData.windInterpolator.vis_interp_phi(lat, lon);
                visualized_values_grid->SetValue(idx_visual, val);
            }
        mapper_grid_main->SetLookupTable(hueLut);
        mapper_grid_main->SetColorModeToMapScalars();
        structuredGrid_main->GetCellData()->SetActiveScalars("visualized_values_grid");

        scalarBar->VisibilityOn();
        scalarBar->SetLookupTable(hueLut);
        scalarBar->SetLabelFormat("%.1f");
        scalarBar_copy1->VisibilityOn();
        scalarBar_copy1->SetLookupTable(hueLut);
        scalarBar_copy1->SetLabelFormat("%.1f");
        hueLut->SetTableRange(0, range);
    }
    else if(VisualizingVariable == VisOpt::wind_from_phi)
    {
        frameData.windInterpolator.prepareVisualizationData(frameData.prms.SimulationTime);
        frameData.windInterpolator.computePotentialApproximation();

        for(int idx_y=0; idx_y<gy; idx_y++)
            for(int idx_x=0; idx_x<gx; idx_x++)
            {
                int idx_visual = idx_x + idx_y*gx;

                float lat = frameData.prms.LatMin + (frameData.prms.LatMax-frameData.prms.LatMin)*(float)idx_y/(float)gy;
                float lon = frameData.prms.LonMin + (frameData.prms.LonMax-frameData.prms.LonMin)*(float)idx_x/(float)gx;

                Eigen::Vector2f wv = frameData.windInterpolator.vis_interp_potential(lat, lon);
                visualized_values_grid->SetValue(idx_visual, wv.norm());
            }
        mapper_grid_main->SetLookupTable(hueLut);
        mapper_grid_main->SetColorModeToMapScalars();
        structuredGrid_main->GetCellData()->SetActiveScalars("visualized_values_grid");

        scalarBar->VisibilityOn();
        scalarBar->SetLookupTable(hueLut);
        scalarBar->SetLabelFormat("%.1f");
        scalarBar_copy1->VisibilityOn();
        scalarBar_copy1->SetLookupTable(hueLut);
        scalarBar_copy1->SetLabelFormat("%.1f");
        hueLut->SetTableRange(0, range);
    }








    int64_t display_date = (int64_t)frameData.prms.SimulationTime + frameData.prms.SimulationStartUnixTime;
    std::time_t unix_time = display_date;
    std::tm* tm_time = std::gmtime(&unix_time);
    // Format the time
    char buffer[100];
    std::strftime(buffer, sizeof(buffer), "%d-%m-%Y %H:%M:%S UTC", tm_time);
    actor_text->SetInput(buffer);
    actor_text_copy1->SetInput(buffer);


}


void VTKVisualization::ChangeVisualizationOption(int option)
{
    VisualizingVariable = (VisOpt)option;
    if(option == VisOpt::P)
    {
        actor_text_title_copy1->SetInput("Hydrostatic pressure");
        actor_text_title->SetInput("Hydrostatic pressure");
    }
    else if(option == VisOpt::Q)
    {
        actor_text_title_copy1->SetInput("Deviatoric stress");
        actor_text_title->SetInput("Deviatoric stress");
    }
    else if(option == VisOpt::Jp_inv)
    {
        actor_text_title_copy1->SetInput("Surf. density");
        actor_text_title->SetInput("Surf. density");
    }
    else if(option == VisOpt::wind_norm)
    {
        actor_text_title_copy1->SetInput("Wind speed");
        actor_text_title->SetInput("Wind speed");
    }
    else
    {
        actor_text_title_copy1->SetInput("-");
        actor_text_title->SetInput("-");
    }

    SynchronizeValues();
}

void VTKVisualization::populateLut(const float lutArray[][3], const int size, vtkNew<vtkLookupTable> &table)
{
    table->SetNumberOfTableValues(size);
    for(int i=0;i<size;i++)
        table->SetTableValue(i,lutArray[i][0],lutArray[i][1],lutArray[i][2]);
    table->SetTableRange(0, size-1);

    table->SetRampToLinear();
}


void VTKVisualization::interpolateLut(const float lutArray[][3], const int size, vtkNew<vtkLookupTable> &table)
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


Eigen::Vector3f VTKVisualization::interpolateColor(const float colorArray[][3], int nColors, float value)
{
    if (nColors < 2)
    {
        spdlog::error("Color array must have at least two colors for interpolation.");
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
    float r_out = (1 - localT) * color1[0] + localT * color2[0];
    float g_out = (1 - localT) * color1[1] + localT * color2[1];
    float b_out = (1 - localT) * color1[2] + localT * color2[2];
    return Eigen::Vector3f(r_out,g_out,b_out);
}
