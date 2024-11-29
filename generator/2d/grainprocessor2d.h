#ifndef GRAINPROCESSOR2D_H
#define GRAINPROCESSOR2D_H

#include <Eigen/Core>
#include "bvhn2d.h"

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
};

#endif // GRAINPROCESSOR2D_H
