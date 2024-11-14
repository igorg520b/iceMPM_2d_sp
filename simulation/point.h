#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath>
#include <utility>
#include <Eigen/Core>

#include "parameters_sim.h"
#include <cuda.h>
#include <cuda_runtime.h>

namespace icy { struct Point; }

struct icy::Point
{
    PointVector2r pos, velocity, vSigma, vSigmaSquared, v_s_hat_tr;
    PointMatrix2r Bp, Fe, U, V; // refer to "The Material Point Method for Simulating Continuum Materials"

    t_PointReal Jp_inv, p_tr, q_tr, Je_tr;
    int grain;
    uint32_t utility_data;
};


#endif // PARTICLE_H
