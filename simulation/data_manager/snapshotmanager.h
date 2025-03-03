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

    void PreparePointsAndSetupGrid(std::string fileName, std::string fileNameModelledArea);   // use PNG image files
    void LoadWindData(std::string fileName);    // netCDF4 data
    void SaveSnapshot(int SimulationStep, double SimulationTime);
    void SaveFrame(int SimulationStep, double SimulationTime);

    void ReadSnapshot(std::string fileName);    // custom HDF5 file

private:
    std::vector<uint8_t> count;   // used for counting points per cell and image generation
    std::vector<uint8_t> rgb;
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
    static void generate_points(int gx, int gy, float points_per_cell, std::vector<std::array<float, 2>> &buffer);


    static constexpr std::array<std::array<float, 3>, 4> colordata_OpenWater = {{
        {0,0,0},
        {0x15/255.,0x1f/255.,0x2f/255.},
        {0x1d/255.,0x29/255.,0x3a/255.},
        {0x28/255.,0x3b/255.,0x52/255.}//,
//        {0x4a/255.,0x4b/255.,0x4b/255.}
    }};

    static constexpr std::array<std::array<float, 3>, 3> colordata_Solid = {{
        {0xac/255.,0xb0/255.,0xb1/255.},
        {0xc4/255.,0xc8/255.,0xcb/255.},
        {0xc8/255.,0xc9/255.,0xca/255.}
    }};

    std::pair<int, float> categorizeColor(const Eigen::Vector3f& rgb);
    static Eigen::Vector3f arrayToEigen(const std::array<float, 3>& arr) {
        return Eigen::Vector3f(arr[0], arr[1], arr[2]);
    }

    static float projectPointOntoCurve(const Eigen::Vector3f& rgb, const std::array<std::array<float, 3>, 3>& curve);
};

#endif // SNAPSHOTWRITER_H
