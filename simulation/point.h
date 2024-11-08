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
    Eigen::Vector2f pos, velocity;
    Eigen::Matrix2f Bp, Fe; // refer to "The Material Point Method for Simulating Continuum Materials"

    float Jp_inv; // track the change in det(Fp)
    short grain;

    float p_tr, q_tr, Je_tr;
    Eigen::Matrix2f U, V;
    Eigen::Vector2f vSigma, vSigmaSquared, v_s_hat_tr;

    uint32_t utility_data;
};


#endif // PARTICLE_H
