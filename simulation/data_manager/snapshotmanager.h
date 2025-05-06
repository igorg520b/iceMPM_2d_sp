#ifndef SNAPSHOTMANAGER_H
#define SNAPSHOTMANAGER_H

#include <array>
#include <vector>
#include <string>
#include <string_view>
#include <filesystem>

#include <H5Cpp.h>
#include <Eigen/Core>

#include "gui/colormap.h"

namespace icy {class SnapshotManager; class Model;}


class icy::SnapshotManager
{
public:
    icy::Model *model;
    std::string SimulationTitle;


    void PrepareGrid(std::string fileNamePNG, std::string fileNameModelledAreaHDF5);
    void PopulatePoints(std::string fileNameModelledAreaHDF5, bool onlyGenerateCache);
    void ReadPointsFromSnapshot(std::string fileNameSnapshotHDF5);
    void SplitIntoPartitionsAndTransferToDevice();


    void SaveSnapshot(int SimulationStep, double SimulationTime);

    void SaveFrame(int SimulationStep, double SimulationTime);



//    void LoadWindData(std::string fileName);    // netCDF4 data
//    void ReadSnapshot(std::string fileName);    // custom HDF5 file

private:

    ColorMap colormap;

    std::vector<uint8_t> count;   // used for counting points per cell and image generation
    std::vector<uint8_t> rgb;
    std::vector<uint8_t> rgb_img_Jpinv, rgb_img, rgb_img_ridges;
    std::vector<float> vis_r, vis_g, vis_b, vis_alpha, vis_Jpinv, vis_P, vis_Q, vis_vx, vis_vy;
    void PrepareFrameArrays(); // invoked from SaveFrame

    constexpr static std::string_view pts_cache_path = "_data/point_cache";

    static constexpr double degreesToRadians(double degrees) { return degrees * M_PI / 180.0; }
    static std::string prepare_file_name(int gx, int gy);
    static bool attempt_to_fill_from_cache(int gx, int gy, std::vector<std::array<float, 2>> &buffer);
    static void generate_and_save(int gx, int gy, float points_per_cell, std::vector<std::array<float, 2>> &buffer);
    static void generate_points(int gx, int gy, float points_per_cell, std::vector<std::array<float, 2>> &buffer);

    static constexpr uint8_t waterColor[3] = {0x15, 0x1f, 0x2f};
    void FillModelledAreaWithBlueColor();
    void SavePointColors();
    void ReadPointColors();

    void SaveImagesJGP(const int frame);
};

#endif // SNAPSHOTWRITER_H
