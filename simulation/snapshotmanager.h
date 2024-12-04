#ifndef SNAPSHOTMANAGER_H
#define SNAPSHOTMANAGER_H

#include <array>
#include <vector>
#include <string>

#include <H5Cpp.h>
#include <Eigen/Core>

namespace icy {class SnapshotManager; class Model;}


class icy::SnapshotManager
{
public:
    icy::Model *model;

    void LoadRawPoints(std::string fileName);

private:
    static constexpr std::array<std::array<float, 3>, 16> colordata = {{
        {0,0,0},                          // water 0
        {0x0c/255.,0x10/255.,0x0f/255.},  // water 1
        {0x0f/255.,0x16/255.,0x1c/255.},  // water 2
        {0x11/255.,0x18/255.,0x20/255.},  // water 3
        {0x17/255.,0x20/255.,0x29/255.},  // water 4
        {0x25/255.,0x39/255.,0x37/255.},  // water 5


        {0x2f/255.,0x42/255.,0x50/255.},  // crushed 6
        {0x3f/255.,0x52/255.,0x56/255.},  // crushed 7
        {0x47/255.,0x5a/255.,0x5e/255.},  // crushed 8
        {0x77/255.,0x9a/255.,0xae/255.},  // crushed 9
        {0xb8/255.,0xca/255.,0xca/255.},  // crushed 10
        {0x81/255.,0xa2/255.,0x99/255.},  // crushed 11
        {0x89/255.,0x76/255.,0xa3/255.},  // crushed 12


        {0xc3/255.,0xc8/255.,0xcb/255.},  // weakened 13

        {0xc7/255.,0xcc/255.,0xcb/255.},  // solid 14
        {1,1,1} // solid 15
    }};

    std::pair<int, float> categorizeColor(const Eigen::Vector3f& rgb);
    static Eigen::Vector3f arrayToEigen(const std::array<float, 3>& arr) {
        return Eigen::Vector3f(arr[0], arr[1], arr[2]);
    }
};

#endif // SNAPSHOTWRITER_H
