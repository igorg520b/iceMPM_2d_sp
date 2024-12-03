#ifndef GRAINPROCESSOR2D_H
#define GRAINPROCESSOR2D_H

#include <Eigen/Core>
#include "bvhn2d.h"
#include <tuple>

struct Triangle
{
    Eigen::Vector2f nds[3];
    int grain;
};


class GrainProcessor2D
{
public:
    std::string outputFileName, meshFileName, landPNGFileName;
    float scale; // for scaling grains
    float block_length;
    float requestedPointsPerCell;
    float cellsize;
    int gridx, gridy;


    void load_png();
    void print_out_parameters();
    void generate_block_and_write();

    std::vector<short> grainID;

    // SOA format (to be written as HDF5)
    std::vector<uint32_t> llGrainID;
    std::vector<float> coordinates[2];

    // mesh with grains
    std::vector<Eigen::Vector2f> vertices2;
    std::vector<std::array<int,4>> elems2;   // 4 nodes + grain id
    std::vector<Triangle> tris2;


private:
    void GenerateBlock();
    void LoadMSH();
    void IdentifyGrains();
    void Write_HDF5();
    static bool PointInsideTriangle(Eigen::Vector2f point, Eigen::Vector2f triangle[3]);

    int nLandNodes, nWaterNodes;
    float volume = -1;
    int imgx, imgy, channels;   // from PNG
    float scale_img;
    std::vector<std::array<float, 2>> buffer;   // result from poisson disk sampler
    std::vector<uint8_t> grid_buffer;
    unsigned char* png_data;

    std::vector<BVHN2D*> leaves;
    BVHN2D root;


    // identification of the point type based on color value
    std::vector<std::tuple<float,float,int>> buffer_categorized;


    static constexpr float colordata[8][3] {
        {0,0,0},// open water
        {0x16/255.,0x1f/255.,0x27/255.},  // open water
        {0x2c/255.,0x3e/255.,0x41/255.},  // crushed
        {0x4e/255.,0x5e/255.,0x60/255.},  // crushed
        {0x6c/255.,0x82/255.,0x8f/255.},  // solid but weakened
        {0x74/255.,0x94/255.,0xa4/255.},  // solid but weakened
        {0xc7/255.,0xcc/255.,0xcb/255.},  // 100% solid
        {1,1,1} // 100% solid
    };

    // Function to calculate the squared distance from a point to a segment
    static float pointToSegmentDistance(const Eigen::Vector3f& p, const Eigen::Vector3f& a, const Eigen::Vector3f& b);

    int categorizeColor(const float color[3]);
};

#endif // GRAINPROCESSOR2D_H
