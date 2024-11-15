#include <cxxopts.hpp>

#include <gmsh.h>
#include <omp.h>
#include <spdlog/spdlog.h>

#include "2d/grainprocessor2d.h"

// -o cov_1m.h5 -m msh_2d/1k_2d.msh -r raster_images/land_600_450b.png

int main(int argc, char *argv[])
{
    gmsh::initialize();
    spdlog::info("testing threads {}", omp_get_max_threads());
#pragma omp parallel
    {     spdlog::info("{}", omp_get_thread_num()); }
    std::cout << std::endl;


    // parse options
    cxxopts::Options options("grain identifier", "Generate raw input file for MPM simulation");

    options.add_options()
        // point generation
        ("o,output", "Output file name", cxxopts::value<std::string>()->default_value("raw_10m.h5"))
        ("m,msh", "Input file with grains", cxxopts::value<std::string>())
        ("r,raster", "Input PNG image file with land map", cxxopts::value<std::string>())

        ("n,pointspercell", "How many points per cell desired", cxxopts::value<float>()->default_value("5.0"))
        ("c,scale", "Scale for grain mapping", cxxopts::value<float>()->default_value("3"))
//        ("r,rescale", "Add a special attribute to the resulting file to rescale later", cxxopts::value<bool>()->default_value("false"))

        // for block
        ("l,length", "Length of the block", cxxopts::value<float>()->default_value("170000"))
        ;

    auto option_parse_result = options.parse(argc, argv);
    GrainProcessor2D gp2d;

    gp2d.outputFileName = option_parse_result["output"].as<std::string>();
    gp2d.meshFileName = option_parse_result["msh"].as<std::string>();
    gp2d.landPNGFileName = option_parse_result["raster"].as<std::string>();

    gp2d.block_length = option_parse_result["length"].as<float>();

    gp2d.scale = option_parse_result["scale"].as<float>();
    gp2d.requestedPointsPerCell = option_parse_result["pointspercell"].as<float>();

    gp2d.generate_block_and_write();
}
