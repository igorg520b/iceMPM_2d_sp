#include "vtk_visualization.h"
#include "parameters_sim.h"
//#include <omp.h>
#include <algorithm>
#include <iostream>
#include <spdlog/spdlog.h>
#include <Eigen/Core>

#include "framedata.h"

VTKVisualization::VTKVisualization(FrameData &fd_) : frameData(fd_)
{
    populateLut(ColorMap::Palette::Pressure, lut_Pressure);
    populateLut(ColorMap::Palette::P2, lut_P2);
    populateLut(ColorMap::Palette::ANSYS, lut_ANSYS);
    populateLut(ColorMap::Palette::Ridges, lut_ridges);

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


void VTKVisualization::ChangeVisualizationOption(int option)
{
    LOGR("VTKVisualization::ChangeVisualizationOption {}", option);
    VisualizingVariable = (VisOpt)option;
    if(option == VisOpt::P)
    {
        actor_text_title->SetInput("Pressure");
    }
    else if(option == VisOpt::Q)
    {
        actor_text_title->SetInput("Dev. stress");
    }
    else if(option == VisOpt::Jp_inv)
    {
        actor_text_title->SetInput("Rel.surf.dens.");
    }
    else if(option == VisOpt::ridges)
    {
        actor_text_title->SetInput("Ridges");
    }
    else if(option == VisOpt::grid_vnorm)
    {
        actor_text_title->SetInput("Vel.mag.");
    }
    else
    {
        actor_text_title->SetInput("-");
    }

    SynchronizeValues();
}



void VTKVisualization::SynchronizeValues()
{
    const SimParams &prms = frameData.ggd.prms;

    const int &width = prms.InitializationImageSizeX;
    const int &height = prms.InitializationImageSizeY;
    const int &ox = prms.ModeledRegionOffsetX;
    const int &oy = prms.ModeledRegionOffsetY;
    const int &gx = prms.GridXTotal;
    const int &gy = prms.GridYTotal;
    const double &h = prms.cellsize;

    if(!width) return;


    renderedImage.assign(frameData.ggd.original_image_colors_rgb.begin(), frameData.ggd.original_image_colors_rgb.end());


    // modify the image


    const double range = std::pow(10,ranges[(int)VisualizingVariable]);
    LOGR("VTKVisualization::SynchronizeValues(), grid [{} x {}], h {}; VV {}; range {:>6.2}", width, height, h,
        VisualizingVariable, range);

#pragma omp parallel for
    for(size_t i=0;i<gx;i++)
        for(size_t j=0;j<gy;j++)
        {
            int status = frameData.ggd.path_indices[(i+ox)+(j+oy)*width];
            if(status == 1000)
            {
                size_t idx2 = j + i*gy;
                std::array<uint8_t, 3> _rgb;
                for(int k=0;k<3;k++) _rgb[k] = frameData.snapshot.rgb[idx2*3+k];


                if(VisualizingVariable == Jp_inv)
                {
                    // visualize Jpinv from the prepared grid / visual array
                    size_t idx2 = j + i*gy;

                    std::array<uint8_t, 3> _rgb;
                    _rgb[0] = frameData.snapshot.rgb[idx2*3+0];
                    _rgb[1] = frameData.snapshot.rgb[idx2*3+1];
                    _rgb[2] = frameData.snapshot.rgb[idx2*3+2];

                    float alpha = frameData.snapshot.mass_mask[idx2] ? 1 : 0;
                    std::array<uint8_t, 3> c = ColorMap::mergeColors(ColorMap::rgb_water, _rgb, alpha);

                    float val = frameData.snapshot.vis_Jpinv[idx2];
                    std::array<uint8_t, 3> c1 = colormap.getColor(ColorMap::Palette::Pressure, 0.5*val/range + 0.5);

                    const float mix_original_color = std::abs(val/range*alpha);

                    std::array<uint8_t, 3> c2 = ColorMap::mergeColors(c, c1, mix_original_color);
                    for(int k=0;k<3;k++) renderedImage[((i+ox) + (j+oy)*width)*3+k] = c2[k];
                }
                if(VisualizingVariable == ridges)
                {
                    // visualize Jpinv from the prepared grid / visual array
                    size_t idx2 = j + i*gy;

                    std::array<uint8_t, 3> _rgb;
                    _rgb[0] = frameData.snapshot.rgb[idx2*3+0];
                    _rgb[1] = frameData.snapshot.rgb[idx2*3+1];
                    _rgb[2] = frameData.snapshot.rgb[idx2*3+2];

                    float alpha = frameData.snapshot.mass_mask[idx2] ? 1 : 0;
                    std::array<uint8_t, 3> c = ColorMap::mergeColors(ColorMap::rgb_water, _rgb, alpha);

                    float val = frameData.snapshot.vis_Jpinv[idx2];
                    std::array<uint8_t, 3> c1 = colormap.getColor(ColorMap::Palette::Ridges, val/range);

                    const float mix_original_color = val > 0 ? alpha * val/range : 0;

                    std::array<uint8_t, 3> c2 = ColorMap::mergeColors(c, c1, mix_original_color);
                    for(int k=0;k<3;k++) renderedImage[((i+ox) + (j+oy)*width)*3+k] = c2[k];
                }

                else if(VisualizingVariable == P)
                {
                    size_t idx2 = j + i*gy;
                    float alpha = frameData.snapshot.mass_mask[idx2] ? 1 : 0;
                    float val = frameData.snapshot.vis_P[idx2];
                    std::array<uint8_t, 3> c = colormap.getColor(ColorMap::Palette::Pressure, 0.5*val/range+0.5);
                    std::array<uint8_t, 3> c2 = ColorMap::mergeColors(ColorMap::rgb_water, c, alpha);
                    for(int k=0;k<3;k++) renderedImage[((i+ox) + (j+oy)*width)*3+k] = c2[k];
                }

                else if(VisualizingVariable == Q)
                {
                    size_t idx2 = j + i*gy;
                    float alpha = frameData.snapshot.mass_mask[idx2] ? 1 : 0;
                    float val = frameData.snapshot.vis_Q[idx2];
                    std::array<uint8_t, 3> c = colormap.getColor(ColorMap::Palette::ANSYS, val/range);
                    std::array<uint8_t, 3> c2 = ColorMap::mergeColors(ColorMap::rgb_water, c, alpha);
                    for(int k=0;k<3;k++) renderedImage[((i+ox) + (j+oy)*width)*3+k] = c2[k];
                }

                else if(VisualizingVariable == grid_vnorm)
                {
                    // visualize Jpinv from the prepared grid / visual array
                    size_t idx2 = j + i*gy;
                    float val_x = frameData.snapshot.vis_vx[idx2];
                    float val_y = frameData.snapshot.vis_vy[idx2];
                    float vnorm = std::sqrt(val_x*val_x + val_y*val_y);
                    std::array<uint8_t, 3> c = colormap.getColor(ColorMap::Palette::ANSYS, vnorm/range);

                    std::array<uint8_t, 3> _rgb;
                    _rgb[0] = frameData.snapshot.rgb[idx2*3+0];
                    _rgb[1] = frameData.snapshot.rgb[idx2*3+1];
                    _rgb[2] = frameData.snapshot.rgb[idx2*3+2];
                    std::array<uint8_t, 3> c2 = ColorMap::mergeColors(_rgb, c, vnorm/range);

                    float alpha = frameData.snapshot.mass_mask[idx2] ? 1 : 0;
                    std::array<uint8_t, 3> c3 = ColorMap::mergeColors(ColorMap::rgb_water, c2, alpha);
                    for(int k=0;k<3;k++) renderedImage[((i+ox) + (j+oy)*width)*3+k] = c3[k];
                }

                else if(VisualizingVariable == colors)
                {
                    // visualize Jpinv from the prepared grid / visual array
                    size_t idx2 = j + i*gy;

                    std::array<uint8_t, 3> _rgb;
                    _rgb[0] = frameData.snapshot.rgb[idx2*3+0];
                    _rgb[1] = frameData.snapshot.rgb[idx2*3+1];
                    _rgb[2] = frameData.snapshot.rgb[idx2*3+2];

                    float alpha = frameData.snapshot.mass_mask[idx2] ? 1 : 0;
                    std::array<uint8_t, 3> c = ColorMap::mergeColors(ColorMap::rgb_water, _rgb, alpha);
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
    raster_imageData->Modified();

    raster_plane->SetOrigin(-h/2, -h/2, -1.0);           // Bottom-left corner
    raster_plane->SetPoint1((width-0.5)*h, -h/2, -1.0);         // Bottom-right (x-axis)
    raster_plane->SetPoint2(-h/2, (height-0.5)*h, -1.0);        // Top-left (y-axis)
    raster_plane->SetNormal(0.0, 0.0, 1.0);           // Normal along z-axis (facing forward)

    raster_mapper->SetInputConnection(raster_plane->GetOutputPort());
    raster_mapper->ScalarVisibilityOff();

    raster_texture->SetInputData(raster_imageData);
    raster_texture->InterpolateOff(); // Smooth texture rendering

    raster_actor->SetMapper(raster_mapper);
    raster_actor->SetTexture(raster_texture);

    raster_mapper->Update();
    raster_texture->Update();

    // update text

    int64_t display_date = (int64_t)frameData.SimulationTime + prms.SimulationStartUnixTime;
    std::time_t unix_time = display_date;
    std::tm* tm_time = std::gmtime(&unix_time);
    // Format the time
    char buffer[100];
    std::strftime(buffer, sizeof(buffer), "%d-%m-%Y %H:%M:%S UTC", tm_time);
    actor_text->SetInput(buffer);

    // scalarbar

    scalarBar->VisibilityOn();
    if(VisualizingVariable == ridges)
    {
        lut_ridges->SetTableRange(0, range);
        scalarBar->SetLookupTable(lut_ridges);
        scalarBar->SetLabelFormat("%.2f");
    }
    else if(VisualizingVariable == Jp_inv)
    {
        lut_Pressure->SetTableRange(-range, range);
        scalarBar->SetLookupTable(lut_Pressure);
        scalarBar->SetLabelFormat("%.1e");
    }
    else if(VisualizingVariable == P)
    {
        lut_Pressure->SetTableRange(-range, range);
        scalarBar->SetLookupTable(lut_Pressure);
        scalarBar->SetLabelFormat("%.1e");
    }
    else if(VisualizingVariable == Q)
    {
        lut_ANSYS->SetTableRange(0, range);
        scalarBar->SetLookupTable(lut_ANSYS);
        scalarBar->SetLabelFormat("%.1e");
    }
    else if(VisualizingVariable == grid_vnorm)
    {
        lut_ANSYS->SetTableRange(0, range);
        scalarBar->SetLookupTable(lut_ANSYS);
        scalarBar->SetLabelFormat("%.2f");
    }
    else
    {
        scalarBar->VisibilityOff();
    }

}



void VTKVisualization::populateLut(ColorMap::Palette palette, vtkNew<vtkLookupTable>& table)
{
    const std::vector<Eigen::Vector3f>& colorTable = ColorMap::getColorTable(palette);
    int size = static_cast<int>(colorTable.size());

    if (size < 2) {
        std::cerr << "Error: Colormap must have at least two colors." << std::endl;
        return;
    }

    const int m = 256;  // Number of colors in the lookup table
    table->SetNumberOfTableValues(m);
    table->Build();

    for (int i = 0; i < m; ++i) {
        float t = static_cast<float>(i) / (m - 1); // Normalize index to [0, 1]

        // Scale t to the range [0, size-1] for interpolation
        float scaledT = t * (size - 1);
        int lowerIdx = static_cast<int>(std::floor(scaledT));
        int upperIdx = static_cast<int>(std::ceil(scaledT));
        float localT = scaledT - lowerIdx; // Fractional part for interpolation

        // Interpolate RGB components
        const Eigen::Vector3f& lowerColor = colorTable[lowerIdx];
        const Eigen::Vector3f& upperColor = colorTable[upperIdx];

        Eigen::Vector3f interpolatedColor = (1.0f - localT) * lowerColor + localT * upperColor;

        // Set interpolated color in the lookup table
        table->SetTableValue(i, interpolatedColor[0], interpolatedColor[1], interpolatedColor[2], 1.0);
    }
}

