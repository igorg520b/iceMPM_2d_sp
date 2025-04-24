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
    colormap.populateLut(ColorMap::Palette::P2, lut_P2);
    colormap.populateLut(ColorMap::Palette::ANSYS, lut_ANSYS);

    // text
    constexpr int fontSize = 20;
    vtkTextProperty* txtprop = actor_text->GetTextProperty();
    txtprop->SetFontFamilyToArial();
    txtprop->BoldOff();
    txtprop->SetFontSize(fontSize);
    txtprop->ShadowOff();
    txtprop->SetColor(0,0,0);

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

    actor_text->SetDisplayPosition(600, 10);

    txtprop = actor_text_title->GetTextProperty();
    txtprop->SetFontFamilyToArial();
    txtprop->BoldOff();
    txtprop->ShadowOff();
    txtprop->SetColor(0,0,0);
    txtprop->SetFontSize(fontSize);

    actor_text_title->SetDisplayPosition(10, 500);
}




void VTKVisualization::SynchronizeValues()
{
    if(!frameData->dataLoaded) return;
    const SimParams &prms = frameData->ggd->prms;

    const int &width = prms.InitializationImageSizeX;
    const int &height = prms.InitializationImageSizeY;
    const int &ox = prms.ModeledRegionOffsetX;
    const int &oy = prms.ModeledRegionOffsetY;
    const int &gx = prms.GridXTotal;
    const int &gy = prms.GridYTotal;
    const double &h = prms.cellsize;



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























    /*
    const SimParams &prms = frameData->ggd->prms;

    const int &_ox = prms.ModeledRegionOffsetX;
    const int &_oy = prms.ModeledRegionOffsetY;

    const int &gx = prms.GridXTotal;
    const int &gy = prms.GridYTotal;

    const double sim_time = frameData->SimulationTime;
    const double &h = prms.cellsize;

    double range = std::pow(10,ranges[VisualizingVariable]);
    double centerVal = 0;
    LOGR("SynchronizeValues() {}; range {}", this->VisualizingVariable, range);

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
        LOGV("rendering Jp_inv");
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
                        float coeff1 = 1;//std::min(count/2.,1.); // how much water surface is obscured
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
        LOGV("rendering colors");
        scalarBar->VisibilityOff();

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
    else if(VisualizingVariable == VisOpt::P)
    {
        for(int i=0;i<gx;i++)
            for(int j=0;j<gy;j++)
            {
                int idx2 = (j + gy*i);
                if(frameData->ggd->grid_status_buffer[idx2])
                {
                    uint8_t count = frameData->count[idx2];
                    float pressure = frameData->vis_P[idx2];
                    std::array<uint8_t, 3> _rgb;
                    for(int k=0;k<3;k++) _rgb[k] = frameData->rgb[idx2*3+k];

                    if(count > 0)
                    {
                        float coeff2 = std::clamp(std::abs(pressure)/range, 0., 1.);
                        coeff2 = pow(coeff2, 2);
                        float coeff1 = 1;//std::min(count/2.,1.); // how much water surface is obscured
                        float val = (pressure)/range + 0.5;
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
        scalarBar->SetLabelFormat("%.1e");
        lut_Pressure->SetTableRange(-range, +range);
    }
    else if(VisualizingVariable == VisOpt::Q)
    {
        for(int i=0;i<gx;i++)
            for(int j=0;j<gy;j++)
            {
                int idx2 = (j + gy*i);
                if(frameData->ggd->grid_status_buffer[idx2])
                {
                    uint8_t count = frameData->count[idx2];
                    float deviatoric_stress = frameData->vis_Q[idx2];
                    std::array<uint8_t, 3> _rgb;
                    for(int k=0;k<3;k++) _rgb[k] = frameData->rgb[idx2*3+k];

                    if(count > 0)
                    {
                        float coeff2 = std::clamp(std::abs(deviatoric_stress)/range, 0., 1.);
                        //coeff2 = pow(coeff2, 0.5);
                        float coeff1 = 1;//std::min(count/2.,1.); // how much water surface is obscured
                        float val = (deviatoric_stress)/range;
                        std::array<uint8_t, 3> c = colormap.getColor(ColorMap::Palette::ANSYS, val);

                        std::array<uint8_t, 3> c2 = ColorMap::mergeColors(_rgb, c, coeff2);
                        std::array<uint8_t, 3> c3 = ColorMap::mergeColors(seaWater, c2, coeff1);

                        unsigned char* pixel = static_cast<unsigned char*>(uniformGrid->GetScalarPointer(i+_ox, j+_oy, 0));
                        for(int k=0;k<3;k++) pixel[k] = c3[k];
                    }
                }
            }
        scalarBar->VisibilityOn();
        scalarBar->SetLookupTable(lut_ANSYS);
        scalarBar->SetLabelFormat("%.1e");
        lut_ANSYS->SetTableRange(0, range);
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

*/
}


void VTKVisualization::ChangeVisualizationOption(int option)
{
    LOGR("VTKVisualization::ChangeVisualizationOption {}", option);
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


