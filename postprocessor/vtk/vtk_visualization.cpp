#include "vtk_visualization.h"
#include "parameters_sim.h"
//#include <omp.h>
#include <algorithm>
#include <iostream>
#include <spdlog/spdlog.h>
#include <Eigen/Core>

#include "framedata.h"

VTKVisualization::VTKVisualization()
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
    if(!frameData->dataLoaded) return;

    const SimParams &prms = frameData->ggd->prms;

    // actor_grid_main
    const int gx = prms.GridXTotal;
    const int gy = prms.GridYTotal;
    const double h = prms.cellsize;

    const int imgx = prms.InitializationImageSizeX;
    const int imgy = prms.InitializationImageSizeY;

    const int ox = prms.ModeledRegionOffsetX;
    const int oy = prms.ModeledRegionOffsetY;

    spdlog::info("img size {}", frameData->ggd->original_image_colors_rgb.size());
    spdlog::info("calculated size {}", imgx*imgy*3);

    uniformGrid->SetDimensions(imgx, imgy, 1);
    uniformGrid->SetSpacing(h, h, 1.0);
    uniformGrid->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    for(int i=0;i<imgx;i++)
        for(int j=0;j<imgy;j++)
        {
            unsigned char* pixel = static_cast<unsigned char*>(uniformGrid->GetScalarPointer(i, j, 0));
            int idx = (i+imgx*j)*3;
            pixel[0] = frameData->ggd->original_image_colors_rgb[idx+0];
            pixel[1] = frameData->ggd->original_image_colors_rgb[idx+1];
            pixel[2] = frameData->ggd->original_image_colors_rgb[idx+2];
        }

    for(int i=0;i<gx;i++)
        for(int j=0;j<gy;j++)
        {
            int idx2 = (j + gy*i);
            if(frameData->ggd->grid_status_buffer[idx2])
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
    if(!frameData->dataLoaded) return;
    const SimParams &prms = frameData->ggd->prms;

    const int &_ox = prms.ModeledRegionOffsetX;
    const int &_oy = prms.ModeledRegionOffsetY;

    const int &gx = prms.GridXTotal;
    const int &gy = prms.GridYTotal;

    const double sim_time = frameData->SimulationTime;
    const double &h = prms.cellsize;
    const double ox = h * prms.ModeledRegionOffsetX;
    const double oy = h * prms.ModeledRegionOffsetY;

    double range = std::pow(10,ranges[VisualizingVariable]);
    double centerVal = 0;
    spdlog::info("SynchronizeValues() {}; range {}", this->VisualizingVariable, range);

    const std::array<uint8_t, 3>& seaWater {0x15, 0x1f, 0x2f};

    for(int i=0;i<gx;i++)
        for(int j=0;j<gy;j++)
        {
            int idx2 = (j + gy*i);
            if(frameData->ggd->grid_status_buffer[idx2])
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
                if(frameData->ggd->grid_status_buffer[idx2])
                {
                    uint8_t count = frameData->count[idx2];
                    float Jp_inv = frameData->vis_Jpinv[idx2];
                    std::array<uint8_t, 3> _rgb;
                    for(int k=0;k<3;k++) _rgb[k] = frameData->rgb[idx2*3+k];

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
                if(frameData->ggd->grid_status_buffer[idx2])
                {
                    uint8_t count = frameData->count[idx2];
                    std::array<uint8_t, 3> _rgb;
                    for(int k=0;k<3;k++) _rgb[k] = frameData->rgb[idx2*3+k];

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

    int64_t display_date = (int64_t)sim_time + frameData->ggd->prms.SimulationStartUnixTime;
    std::time_t unix_time = display_date;
    std::tm* tm_time = std::gmtime(&unix_time);
    // Format the time
    char buffer[100];
    std::strftime(buffer, sizeof(buffer), "%d-%m-%Y %H:%M:%S UTC", tm_time);
    actor_text->SetInput(buffer);
}


void VTKVisualization::ChangeVisualizationOption(int option)
{
    spdlog::info("VTKVisualization::ChangeVisualizationOption {}", option);
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

    SynchronizeTopology();
}


