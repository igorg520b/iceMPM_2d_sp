#ifndef SNAPSHOTMANAGER_H
#define SNAPSHOTMANAGER_H

#include <array>
#include <vector>
#include <string>
#include <string_view>
#include <filesystem>


#include <H5Cpp.h>
#include <Eigen/Core>

namespace icy {class SnapshotManager; class Model;}


class icy::SnapshotManager
{
public:
    icy::Model *model;

    void PreparePointsAndSetupGrid(std::string fileName);   // use PNG image file
    void LoadWindData(std::string fileName);    // netCDF4 data
    void SaveSnapshot(int SimulationStep, double SimulationTime);
    void SaveFrame(int SimulationStep, double SimulationTime);

private:
    std::vector<uint8_t> tmp;   // used for counting points per cell and image generation
    std::vector<float> vis_r, vis_g, vis_b, vis_alpha, vis_Jpinv, vis_P, vis_Q, vis_vx, vis_vy;

    constexpr static std::string_view pts_cache_path = "_data/point_cache";
    constexpr static std::string_view snapshot_path = "_data/snapshots";
    constexpr static std::string_view frame_path = "_data/frames";
    constexpr static std::string_view image_path = "_data/images";

    static constexpr double degreesToRadians(double degrees) { return degrees * M_PI / 180.0; }
    static double haversineDistance(double lat, double lon1, double lon2);
    static std::string prepare_file_name(int gx, int gy);
    static bool attempt_to_fill_from_cache(int gx, int gy, std::vector<std::array<float, 2>> &buffer);
    static void generate_and_save(int gx, int gy, float points_per_cell, std::vector<std::array<float, 2>> &buffer);

    void load_png(std::string pngFileName, unsigned char* &png_data);


    static constexpr std::array<std::array<float, 3>, 6> colordata_OpenWater = {{
        {0,0,0},
        {0x0c/255.,0x10/255.,0x0f/255.},
        {0x0f/255.,0x16/255.,0x1c/255.},
        {0x11/255.,0x18/255.,0x20/255.},
        {0x17/255.,0x20/255.,0x29/255.},
        {0x25/255.,0x39/255.,0x37/255.},
    }};

    static constexpr std::array<std::array<float, 3>, 8> colordata_Crushed = {{
        {0x2f/255.,0x42/255.,0x50/255.}, // 1
        {0x3f/255.,0x52/255.,0x56/255.}, // 2
        {0x47/255.,0x5a/255.,0x5e/255.}, // 3
        {0x77/255.,0x9a/255.,0xae/255.}, // 4
        {0x89/255.,0x76/255.,0xa3/255.}, // 7
        {0x81/255.,0xa2/255.,0x99/255.}, // 6
        {0xc3/255.,0xc8/255.,0xcb/255.}, // 8
        {0xb8/255.,0xca/255.,0xca/255.}, // 5
    }};

    std::pair<int, float> categorizeColor(const Eigen::Vector3f& rgb);
    static Eigen::Vector3f arrayToEigen(const std::array<float, 3>& arr) {
        return Eigen::Vector3f(arr[0], arr[1], arr[2]);
    }

    static float projectPointOntoCurve(const Eigen::Vector3f& rgb, const std::array<std::array<float, 3>, 8>& curve);
};

#endif // SNAPSHOTWRITER_H
