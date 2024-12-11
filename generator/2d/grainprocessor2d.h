#ifndef GRAINPROCESSOR2D_H
#define GRAINPROCESSOR2D_H

#include <Eigen/Core>
#include <tuple>



class GrainProcessor2D
{
public:
    std::string outputFileName, landPNGFileName;
    float requestedPointsPerCell;
    float cellsize;
    int gridx, gridy;


    void load_png();
    void print_out_parameters();
    void generate_block_and_write();

//    std::vector<short> grainID;

    // SOA format (to be written as HDF5)
    std::vector<uint32_t> llGrainID;
    std::vector<float> coordinates[2];
    std::vector<float> rgb[3];

private:
    void GenerateBlock();
    void IdentifyGrains();
    void Write_HDF5();

    int nLandNodes, nWaterNodes;
    int imgx, imgy, channels;   // from PNG
    float scale_img;
    std::vector<uint8_t> grid_buffer;
    unsigned char* png_data;


};

#endif // GRAINPROCESSOR2D_H
