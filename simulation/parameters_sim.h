#ifndef P_SIM_H
#define P_SIM_H

#include <iostream>
#include <string>
#include <filesystem>
#include <fstream>

#include <Eigen/Core>
#include "rapidjson/reader.h"
#include "rapidjson/document.h"
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

// variables related to the formulation of the model

typedef double t_GridReal;      // data type for grid data
typedef double t_PointReal;     // data type to store point data
//typedef float t_GridReal;      // data type for grid data
//typedef float t_PointReal;     // data type to store point data

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
    constexpr static int dim = 2;
    constexpr static int nGridArrays = 3; // mass, px, py

    // index of the corresponding array in SoA
    constexpr static size_t idx_utility_data = 0;
    constexpr static size_t integer_cell_idx = idx_utility_data + 1;
    constexpr static size_t idx_P = integer_cell_idx + 1;
    constexpr static size_t idx_Q = idx_P + 1;
    constexpr static size_t idx_Jp_inv = idx_Q + 1;
    constexpr static size_t idx_Qp = idx_Jp_inv + 1;    // accumulated plastic shear
    constexpr static size_t posx = idx_Qp + 1;
    constexpr static size_t velx = posx + 2;
    constexpr static size_t Fe00 = velx + 2;
    constexpr static size_t Bp00 = Fe00 + 4;
    constexpr static size_t nPtsArrays = Bp00 + 4;

    int tpb_P2G, tpb_Upd, tpb_G2P;  // threads per block for each operation

    int nPtsInitial;
    int GridXTotal, GridY;

    t_PointReal InitialTimeStep, SimulationEndTime;
    int AnimationFramesRequested, UpdateEveryNthStep; // run N steps without update
    int SimulationStep;
    t_PointReal SimulationTime;

    // wind
    t_PointReal windDragCoeff_airDensity, waterDrag_waterDensity;
    t_PointReal GlenA;

    // material properties
    t_PointReal Gravity, Density, PoissonsRatio, YoungsModulus;
    t_PointReal lambda, mu, kappa; // Lame

    t_PointReal IceCompressiveStrength, IceTensileStrength, IceShearStrength, IceTensileStrength2;
    t_PointReal NACC_beta, NACC_M, NACC_Msq;     // these are all computed

    t_PointReal DP_tan_phi, DP_threshold_p;

    // indentation params
    t_PointReal IndDiameter, IndRSq, IndVelocity, IndDepth;

    t_PointReal cellsize, cellsize_inv, Dp_inv;
    t_PointReal xmin, xmax, ymin, ymax;            // bounding box of the material
    int nxmin, nxmax, nymin, nymax;         // same, but nuber of grid cells

    t_PointReal ParticleVolume, ParticleMass, ParticleViewSize;

    t_PointReal Volume;  // total volume (area) of the object
    int SetupType;  // 0 - ice block horizontal indentation; 1 - cone uniaxial compression
    t_PointReal GrainVariability;

    // computed parameters/properties
    t_PointReal dt_vol_Dpinv, dt_Gravity, vmax, vmax_squared;

    void Reset();
    std::string ParseFile(std::string fileName);

    void ComputeLame();
    void ComputeCamClayParams2();
    void ComputeHelperVariables();
    void ComputeIntegerBlockCoords();
    float PointsPerCell() {return nPtsInitial/(Volume/(cellsize*cellsize));}
    int AnimationFrameNumber() { return SimulationStep / UpdateEveryNthStep;}
};

#endif
