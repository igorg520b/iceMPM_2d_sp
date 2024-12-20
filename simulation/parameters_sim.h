#ifndef P_SIM_H
#define P_SIM_H

#include <iostream>
#include <string>
#include <filesystem>
#include <fstream>
#include <map>

#include <Eigen/Core>
#include "rapidjson/reader.h"
#include "rapidjson/document.h"
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include <H5Cpp.h>
// variables related to the formulation of the model

typedef double real;
//typedef float real;
typedef real t_GridReal;      // data type for grid data
typedef real t_PointReal;     // data type to store point data

typedef Eigen::Matrix<t_GridReal, 2, 1> GridVector2r;
typedef Eigen::Matrix<t_GridReal, 2, 2> GridMatrix2r;
typedef Eigen::Matrix<t_PointReal, 2, 1> PointVector2r;
typedef Eigen::Matrix<t_PointReal, 2, 2> PointMatrix2r;
typedef Eigen::Array<t_PointReal, 2, 1> PointArray2r;

struct SimParams
{
public:
    constexpr static float disabled_pts_proportion_threshold = 0.05; // when exceeded, disabled points are removed
    constexpr static t_PointReal pi = 3.14159265358979323846;
    constexpr static double Earth_Radius = 6371000.0;
    constexpr static float MPM_points_per_cells = 5.0;    // approximate average value

    constexpr static int dim = 2;
    constexpr static int nGridArrays = 3; // mass, px, py
    constexpr static int snapshotFrequency = 25;

    // index of the corresponding array in SoA
    constexpr static size_t idx_utility_data = 0;
    constexpr static size_t integer_cell_idx = idx_utility_data + 1;
    constexpr static size_t integer_point_idx = integer_cell_idx + 1;
    constexpr static size_t idx_P = integer_point_idx + 1;
    constexpr static size_t idx_Q = idx_P + 1;
    constexpr static size_t idx_Jp_inv = idx_Q + 1;
    constexpr static size_t idx_Qp = idx_Jp_inv + 1;    // accumulated plastic shear
    constexpr static size_t posx = idx_Qp + 1;
    constexpr static size_t velx = posx + 2;
    constexpr static size_t Fe00 = velx + 2;
    constexpr static size_t Bp00 = Fe00 + 4;
    constexpr static size_t idx_initial_strength = Bp00 + 4;
    constexpr static size_t nPtsArrays = idx_initial_strength + 1;

    int tpb_P2G, tpb_Upd, tpb_G2P;  // threads per block for each operation

    int nPtsInitial;
    int64_t SimulationStartUnixTime;
    int GridXTotal, GridY;
    double LatMin, LatMax, LonMin, LonMax; // coordinates corresponding to the modelled space
    double gridLatMin, gridLonMin;
    double DimensionHorizontal, DimensionVertical;

    double InitialTimeStep, SimulationEndTime;
    int AnimationFramesRequested; // run N steps without update
    int SimulationStep;
    double SimulationTime;

    // wind
    double windDragCoeff_airDensity;

    // material properties
    double SurfaceDensity, PoissonsRatio, YoungsModulus;
    double IceCompressiveStrength, IceTensileStrength, IceShearStrength, IceTensileStrength2;
    double DP_phi, DP_threshold_p;
    double cellsize;
    double ParticleVolume, ParticleViewSize;

    double GrainVariability;

    // computed parameters/properties
    double dt_vol_Dpinv, vmax, vmax_squared;
    double lambda, mu, kappa; // Lame
    double ParticleMass;
    double cellsize_inv, Dp_inv;
    int UpdateEveryNthStep;

    void Reset();
    std::map<std::string,std::string> ParseFile(std::string fileName);  // return additional filenames to load

    void ComputeLame();
    void ComputeHelperVariables();
    int AnimationFrameNumber() { return SimulationStep / UpdateEveryNthStep;}

    void SaveParametersAsHDF5Attributes(H5::DataSet &dataset);
    void ReadParametersFromHDF5Attributes(H5::DataSet &dataset);
};

#endif
