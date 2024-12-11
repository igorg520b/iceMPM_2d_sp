#include "grainprocessor2d.h"
#include "poisson_disk_sampling.h"

#include <vector>
#include <unordered_map>
#include <filesystem>
#include <algorithm>
#include <iostream>
#include <algorithm>

#include <H5Cpp.h>
#include <gmsh.h>
#include <spdlog/spdlog.h>

#include <Eigen/Dense>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"



void GrainProcessor2D::generate_block_and_write()
{
    load_png();
    print_out_parameters();
    GenerateBlock();
    Write_HDF5();
    stbi_image_free(png_data);
}




void GrainProcessor2D::print_out_parameters()
{
    spdlog::info("requestedPointsPerCell {}", requestedPointsPerCell);
    spdlog::info("grid dimensions {} x {}", gridx, gridy);
    spdlog::info("proportion of landmass {}",(float)nLandNodes/(gridx*gridy));
    spdlog::info("expected points {}", (gridx*gridy-nLandNodes)*requestedPointsPerCell);
    spdlog::info("output file: {}", outputFileName);
}


void GrainProcessor2D::load_png()
{
    spdlog::info("load_png");
    const char* filename = landPNGFileName.c_str();
    png_data = stbi_load(filename, &imgx, &imgy, &channels, 3);  // Request 1 channel (grayscale)

    if (!png_data)
    {
        std::cerr << "Failed to load image: " << filename << std::endl;
        throw std::runtime_error("png not loaded");
    }

    if(gridx == -1) gridx = imgx;

    scale_img = (float)imgx / gridx;
    gridy = (int)(imgy/scale_img);
    grid_buffer.resize(gridx*gridy);
    spdlog::info("image ({} x {}) loaded; {} channels", imgx, imgy, channels);
    spdlog::info("grid size {} x {}", gridx, gridy);
    spdlog::info("grid_buffer size {}",grid_buffer.size());

    nLandNodes = 0;
    nWaterNodes = 0;
    for (int y = 0; y < gridy; ++y)
    {
        for (int x = 0; x < gridx; ++x)
        {
            int px = (int)(x * scale_img);
            int py = (int)(y * scale_img);
            px = std::min(px, imgx-1);
            py = std::min(py, imgy-1);

            int index_png = py * imgx + px;
            unsigned char &r = png_data[index_png*3 + 0];   // red channel
            unsigned char &g = png_data[index_png*3 + 1];   // red channel
            unsigned char &b = png_data[index_png*3 + 2];   // red channel
            bool is_land = (r > 128 && g < 128 & b < 128);
            if(!is_land) nWaterNodes++;
            else nLandNodes++;

            grid_buffer[(gridy-y-1) + x*gridy] = is_land ? 1 : 0;
        }
    }
}




void GrainProcessor2D::GenerateBlock()
{
    constexpr float magic_constant = 0.61;

    float h_img = 1.f/(float(imgx-1));

    float h = 1.f/(float)(gridx-1);
    float dx = 1.f;
    float dy = 1.f*(gridy-1)/(float)(gridx-1);
    float nPts = requestedPointsPerCell*gridx*gridy;
    const float kRadius = sqrt(magic_constant*dx*dy/nPts);

    const std::array<float, 2>kXMin {2*h, 2*h};
    const std::array<float, 2>kXMax {dx-2*h, dy-2*h};

    std::vector<std::array<float, 2>> buffer = thinks::PoissonDiskSampling(kRadius, kXMin, kXMax);
    spdlog::info("generating rectangle [{},{}] - [{},{}]",kXMin[0],kXMin[1],kXMax[0],kXMax[1]);

    // copy to buffer_categorized and categorize
    size_t n = buffer.size();
    coordinates[0].resize(n);
    coordinates[1].resize(n);
    rgb[0].resize(n);
    rgb[1].resize(n);
    rgb[2].resize(n);
    llGrainID.resize(n);

    float max_x = 0, max_y = 0;
    for(int idx=0;idx<buffer.size();idx++)
    {
        float x = buffer[idx][0];
        float y = buffer[idx][1];

        int i = (int)(x/h + 0.5f);
        int j = (int)(y/h + 0.5f);
        if(i<0 || j< 0 || i>=gridx || j>=gridy)
        {
            spdlog::critical("error pt ({},{}); grid cell ({},{})",x,y,i,j);
            throw std::runtime_error("particle is out of grid bounds");
        }

        int px = (int)(x/h_img + 0.5f);
        int py = (int)(y/h_img + 0.5f);
        int index_png = (imgy - py -1) * imgx + px;
        rgb[0][idx] = (float)png_data[index_png*3 + 0]/255.;
        rgb[1][idx] = (float)png_data[index_png*3 + 1]/255.;
        rgb[2][idx] = (float)png_data[index_png*3 + 2]/255.;
        coordinates[0][idx] = x;
        coordinates[1][idx] = y;
        max_x = std::max(max_x, coordinates[0][idx]);
        max_y = std::max(max_y, coordinates[1][idx]);
    }

    spdlog::info("final extent of the block: x {}; y {}", max_x, max_y);
    spdlog::info("points per cell = {}", (float)n/(gridx*gridy));

}


void GrainProcessor2D::Write_HDF5()
{
    spdlog::info("Write_HDF5()");
    H5::H5File file(outputFileName, H5F_ACC_TRUNC);

    const size_t n = coordinates[0].size();
    hsize_t dims_pts[1] = {n};
    H5::DataSpace dataspace_pts(1, dims_pts);

    hsize_t g_chunk_dims[1] = {1024*256};
    if(g_chunk_dims[0] > n) g_chunk_dims[0] = std::max(n/5,(size_t)1);
    spdlog::info("chunk {}; dims_pts {}", g_chunk_dims[0], dims_pts[0]);
    H5::DSetCreatPropList proplist1;
    proplist1.setChunk(1, g_chunk_dims);
    proplist1.setDeflate(7);

    H5::DataSet dataset_x = file.createDataSet("x", H5::PredType::NATIVE_FLOAT, dataspace_pts, proplist1);
    dataset_x.write(coordinates[0].data(), H5::PredType::NATIVE_FLOAT);

    H5::DataSet dataset_y = file.createDataSet("y", H5::PredType::NATIVE_FLOAT, dataspace_pts, proplist1);
    dataset_y.write(coordinates[1].data(), H5::PredType::NATIVE_FLOAT);

    H5::DataSet dataset_r = file.createDataSet("r", H5::PredType::NATIVE_FLOAT, dataspace_pts, proplist1);
    dataset_r.write(rgb[0].data(), H5::PredType::NATIVE_FLOAT);

    H5::DataSet dataset_g = file.createDataSet("g", H5::PredType::NATIVE_FLOAT, dataspace_pts, proplist1);
    dataset_g.write(rgb[1].data(), H5::PredType::NATIVE_FLOAT);

    H5::DataSet dataset_b = file.createDataSet("b", H5::PredType::NATIVE_FLOAT, dataspace_pts, proplist1);
    dataset_b.write(rgb[2].data(), H5::PredType::NATIVE_FLOAT);


    hsize_t att_dim = 1;
    H5::DataSpace att_dspace(1, &att_dim);

    // write grid info
    hsize_t dims_grid[1] = {grid_buffer.size()};

    H5::DataSpace dataspace_grid(1, dims_grid);
    hsize_t g_chunk_grid[1] = {1024*128};
    if(g_chunk_grid[0] > grid_buffer.size()) g_chunk_grid[0] = std::max(grid_buffer.size()/3,(size_t)1);
    H5::DSetCreatPropList proplist2;

    proplist2.setChunk(1, g_chunk_grid);
    proplist2.setDeflate(7);
    H5::DataSet dataset_grid = file.createDataSet("GridLand", H5::PredType::NATIVE_UINT8, dataspace_grid, proplist2);
    dataset_grid.write(grid_buffer.data(), H5::PredType::NATIVE_UINT8);

    dataset_grid.createAttribute("GridX", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &gridx);
    dataset_grid.createAttribute("GridY", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &gridy);

    dataset_grid.createAttribute("nLandNodes", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &nLandNodes);
    dataset_grid.createAttribute("nWaterNodes", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &nWaterNodes);

    file.close();
}

