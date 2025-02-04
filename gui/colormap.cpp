#include "colormap.h"


void ColorMap::populateLut(Palette palette, vtkNew<vtkLookupTable> &table)
{
    const auto& colors = colormaps[static_cast<size_t>(palette)];
    int size = static_cast<int>(colors.size());

    table->SetNumberOfTableValues(size);

    for (int i = 0; i < size; i++) {
        Eigen::Vector3f color = colors[i];
        table->SetTableValue(i, static_cast<double>(color.x()),
                             static_cast<double>(color.y()),
                             static_cast<double>(color.z()));
    }

    table->SetTableRange(0.0, 1.0);  // Normalized range
    table->SetRampToLinear();
}


Eigen::Vector3f ColorMap::interpolateColor(Palette palette, float value)
{
    const auto& colors = colormaps[static_cast<size_t>(palette)];
    if (colors.empty()) return {0, 0, 0}; // Fallback in case of empty palette

    value = std::clamp(value, 0.0f, 1.0f);
    float scaled = value * (colors.size() - 1);
    int idx = static_cast<int>(scaled);
    float t = scaled - idx;

    if (idx >= colors.size() - 1) return colors.back();
    return (1.0f - t) * colors[idx] + t * colors[idx + 1];
}

// Get interpolated color as uint8_t[3] (values in range [0,255])
std::array<uint8_t, 3> ColorMap::getColor(Palette palette, float value)
{
    Eigen::Vector3f color = interpolateColor(palette, value) * 255.0f;
    return { static_cast<uint8_t>(color.x()), static_cast<uint8_t>(color.y()), static_cast<uint8_t>(color.z()) };
}

// Get full color table for a given palette
const std::vector<Eigen::Vector3f>& ColorMap::getColorTable(Palette palette)
{
    return colormaps[static_cast<size_t>(palette)];
}

std::array<uint8_t, 3> ColorMap::mergeColors(uint32_t rgb, const std::array<uint8_t, 3>& colorArray, float alpha)
{
    // Ensure alpha is clamped between 0 and 1
    alpha = std::max(0.0f, std::min(1.0f, alpha));

    // Extract RGB components from the uint32_t color
    uint8_t r1 = (rgb >> 16) & 0xFF;
    uint8_t g1 = (rgb >> 8) & 0xFF;
    uint8_t b1 = rgb & 0xFF;

    // Get the second color components
    uint8_t r2 = colorArray[0];
    uint8_t g2 = colorArray[1];
    uint8_t b2 = colorArray[2];

    // Perform linear interpolation for each channel
    uint8_t r = static_cast<uint8_t>(r1 * (1 - alpha) + r2 * alpha);
    uint8_t g = static_cast<uint8_t>(g1 * (1 - alpha) + g2 * alpha);
    uint8_t b = static_cast<uint8_t>(b1 * (1 - alpha) + b2 * alpha);

    // Return the result as an std::array<uint8_t, 3>
    return {r, g, b};
}


const std::array<std::vector<Eigen::Vector3f>, static_cast<size_t>(ColorMap::Palette::COUNT)>
    ColorMap::colormaps = {{
    // SpecialJ
    {
        {0.342992, 0.650614, 0.772702},
        {0.688376, 0.931066, 0.963615},
        {0.836049, 0.882901, 0.85903},
        {0xee/255.,0x6e/255.,0xba/255.},
        {0x94/255.,0x00/255.,0x58/255.}
    },

                            {{0xed/255.,0xdf/255.,0xd6/255.},
                            {0xc6/255.,0xb8/255.,0xaf/255.},
                            {0x9c/255.,0x8b/255.,0x7b/255.},
                            {0x7a/255.,0x63/255.,0x55/255.},
                            {0x52/255.,0x40/255.,0x32/255.},
                             {0x2d/255.,0x1e/255.,0x17/255.}},

    // VIRIDIS colormap (9 colors)
    {
        {0.267, 0.004, 0.329}, {0.283, 0.141, 0.458}, {0.254, 0.265, 0.530},
        {0.207, 0.372, 0.553}, {0.164, 0.471, 0.558}, {0.133, 0.567, 0.553},
        {0.122, 0.659, 0.537}, {0.220, 0.741, 0.502}, {0.477, 0.821, 0.318}
    },
    // PLASMA colormap (5 colors)
    {
        {0.050, 0.029, 0.528}, {0.204, 0.027, 0.662}, {0.390, 0.072, 0.705},
        {0.586, 0.165, 0.675}, {0.768, 0.247, 0.580}
    }
}};


