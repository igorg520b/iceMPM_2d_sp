#include <cuda.h>
#include <cuda_runtime.h>
#include "point.h"



using namespace Eigen;

constexpr float d = 2; // dimensions
constexpr float coeff1 = 1.4142135623730950; // sqrt((6-d)/2.);
constexpr float coeff1_inv = 0.7071067811865475;
constexpr uint32_t status_crushed = 0x10000;
constexpr uint32_t status_disabled = 0x20000;

__device__ uint8_t gpu_error_indicator;
__constant__ icy::SimParams gprms;


__device__ Matrix2f KirchhoffStress_Wolper(const Matrix2f &F)
{
    const float &kappa = gprms.kappa;
    const float &mu = gprms.mu;

    // Kirchhoff stress as per Wolper (2019)
    float Je = F.determinant();
    Matrix2f b = F*F.transpose();
    Matrix2f PFt = mu*(1/Je)*dev(b) + kappa*0.5f*(Je*Je-1.f)*Matrix2f::Identity();
    return PFt;
}


__global__ void partition_kernel_p2g(const int gridX, const int gridX_offset, const int pitch_grid,
                                     const int count_pts, const int pitch_pts,
                                     const float *buffer_pts, float *buffer_grid)
{
    int pt_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(pt_idx >= count_pts) return;

    const uint32_t* ptr = reinterpret_cast<const uint32_t*>(&buffer_pts[pitch_pts*icy::SimParams::idx_utility_data]);
    uint32_t utility_data = ptr[pt_idx];
    if(utility_data & status_disabled) return; // point is disabled

    const float &h = gprms.cellsize;
    const float &h_inv = gprms.cellsize_inv;
    const float particle_mass = gprms.ParticleMass;

    const int &gridY = gprms.GridY;

    // pull point data from SOA
    Vector2f pos, velocity;
    Matrix2f Bp, Fe;

    for(int i=0; i<icy::SimParams::dim; i++)
    {
        pos[i] = buffer_pts[pt_idx + pitch_pts*(icy::SimParams::posx+i)];
        velocity[i] = buffer_pts[pt_idx + pitch_pts*(icy::SimParams::velx+i)];
        for(int j=0; j<icy::SimParams::dim; j++)
        {
            Fe(i,j) = buffer_pts[pt_idx + pitch_pts*(icy::SimParams::Fe00 + i*icy::SimParams::dim + j)];
            Bp(i,j) = buffer_pts[pt_idx + pitch_pts*(icy::SimParams::Bp00 + i*icy::SimParams::dim + j)];
        }
    }

    Matrix2f PFt;
    PFt = KirchhoffStress_Wolper(Fe);

    Matrix2f subterm2 = particle_mass*Bp - (gprms.dt_vol_Dpinv)*PFt;

    Eigen::Vector2i base_coord_i = (pos*h_inv - Vector2f::Constant(0.5f)).cast<int>(); // coords of base grid node for point
    Vector2f base_coord = base_coord_i.cast<float>();
    Vector2f fx = pos*h_inv - base_coord;

    if(base_coord_i.x() - gridX_offset < 0) gpu_error_indicator = 70;
    if(base_coord_i.y()<0) gpu_error_indicator = 71;
    if(base_coord_i.x() - gridX_offset>(gridX-3)) gpu_error_indicator = 72;
    if(base_coord_i.y()>gridY-3) gpu_error_indicator = 73;

    // optimized method of computing the quadratic (!) weight function (no conditional operators)
    Array2f arr_v0 = 1.5f-fx.array();
    Array2f arr_v1 = fx.array() - 1.0f;
    Array2f arr_v2 = fx.array() - 0.5f;
    Array2f ww[3] = {0.5f*arr_v0*arr_v0, 0.75f-arr_v1*arr_v1, 0.5f*arr_v2*arr_v2};

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            float Wip = ww[i][0]*ww[j][1];
            Vector2f dpos((i-fx[0])*h, (j-fx[1])*h);
            Vector2f incV = Wip*(velocity*particle_mass + subterm2*dpos);
            float incM = Wip*particle_mass;

            // the x-index of the cell takes into accout the partition's offset of the gird fragment
            int i2 = i+base_coord_i[0] - gridX_offset;
            int j2 = j+base_coord_i[1];
            int idx_gridnode = j2 + i2*gridY;
            // Udpate mass, velocity and force
            atomicAdd(&buffer_grid[0*pitch_grid + idx_gridnode], incM);
            atomicAdd(&buffer_grid[1*pitch_grid + idx_gridnode], incV[0]);
            atomicAdd(&buffer_grid[2*pitch_grid + idx_gridnode], incV[1]);
        }
}




__global__ void partition_kernel_update_nodes(const Eigen::Vector2f indCenter,
                                              const int nNodes, const int gridX_offset, const int pitch_grid,
                                              float *buffer_grid, float *indenter_force_accumulator,
                                              float simulation_time)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= nNodes) return;

    const int &gridY = gprms.GridY;

    float mass = buffer_grid[idx];
    float vx = buffer_grid[1*pitch_grid + idx];
    float vy = buffer_grid[2*pitch_grid + idx];
    if(mass == 0) return;

    const float &indRsq = gprms.IndRSq;
    const float &dt = gprms.InitialTimeStep;
    const float &ind_velocity = gprms.IndVelocity;
    const float &cellsize = gprms.cellsize;
    const float &vmax = gprms.vmax;
    const float &vmax_squared = gprms.vmax_squared;
    const int &gridXTotal = gprms.GridXTotal;

    const Vector2f vco(ind_velocity,0);  // velocity of the collision object (indenter)

    Vector2i gi(idx/gridY+gridX_offset, idx%gridY);   // integer x-y index of the grid node
    Vector2f gnpos = gi.cast<float>()*cellsize;    // position of the grid node in the whole grid

    Vector2f velocity(vx, vy);
    velocity /= mass;
    velocity[1] -= gprms.dt_Gravity;



    if(velocity.squaredNorm() > vmax_squared) velocity = velocity.normalized()*vmax;


    // indenter
    Vector2f n = gnpos - indCenter;
    if(n.squaredNorm() < indRsq)
    {
        // grid node is inside the indenter
        Vector2f vrel = velocity - vco;
        n.normalize();
        float vn = vrel.dot(n);   // normal component of the velocity
        if(vn < 0)
        {
            Vector2f vt = vrel - n*vn;   // tangential portion of relative velocity
            Vector2f prev_velocity = velocity;
            velocity = vco + vt;

            // force on the indenter
            Vector2f force = (prev_velocity-velocity)*mass/dt;
            float angle = atan2f((float)n[0],(float)n[1]);
            angle += icy::SimParams::pi;
            angle *= gprms.n_indenter_subdivisions/ (2*icy::SimParams::pi);
            int index = min(max((int)angle, 0), gprms.n_indenter_subdivisions-1);
            atomicAdd(&indenter_force_accumulator[0+2*index], force[0]);
            atomicAdd(&indenter_force_accumulator[1+2*index], force[1]);
        }
    }

    // attached bottom layer
    if(gi.y() <= 2) velocity.setZero();
//    if(gi.y() <= 2 && velocity[1]<0) velocity[1] = 0;
    if(gi.y() >= gridY-3 && velocity[1]>0) velocity[1] = 0;
    if(gi.x() <= 2 && velocity[0]<0) velocity[0] = 0;
    else if(gi.x() >= gridXTotal-3 && velocity[0]>0) velocity[0] = 0;

    // side boundary conditions would go here

    // write the updated grid velocity back to memory
    buffer_grid[1*pitch_grid + idx] = velocity[0];
    buffer_grid[2*pitch_grid + idx] = velocity[1];
}



__global__ void partition_kernel_g2p(const bool recordPQ, const bool enablePointTransfer,
                                     const int gridX, const int gridX_offset, const int pitch_grid,
                                     const int count_pts, const int pitch_pts,
                                     float *buffer_pts, const float *buffer_grid)
{
    const int pt_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(pt_idx >= count_pts) return;

    // skip if a point is disabled
    icy::Point p;
    uint32_t* ptr = reinterpret_cast<uint32_t*>(&buffer_pts[pt_idx + pitch_pts*icy::SimParams::idx_utility_data]);
    p.utility_data = *ptr;
    if(p.utility_data & status_disabled) return; // point is disabled

    const float &h_inv = gprms.cellsize_inv;
    const float &dt = gprms.InitialTimeStep;
    const int &gridY = gprms.GridY;
    const float &mu = gprms.mu;
    const float &kappa = gprms.kappa;

    // pull point data from SOA
    for(int i=0; i<icy::SimParams::dim; i++)
    {
        p.pos[i] = buffer_pts[pt_idx + pitch_pts*(icy::SimParams::posx+i)];
        for(int j=0; j<icy::SimParams::dim; j++)
        {
            p.Fe(i,j) = buffer_pts[pt_idx + pitch_pts*(icy::SimParams::Fe00 + i*icy::SimParams::dim + j)];
        }
    }
    p.Jp_inv = buffer_pts[pt_idx + pitch_pts*icy::SimParams::idx_Jp_inv];
    p.grain = (short)p.utility_data;

    // coords of base grid node for point
    Eigen::Vector2i base_coord_i = (p.pos*h_inv - Vector2f::Constant(0.5f)).cast<int>();
    Vector2f base_coord = base_coord_i.cast<float>();
    Vector2f fx = p.pos*h_inv - base_coord;

    // optimized method of computing the quadratic weight function without conditional operators
    Array2f arr_v0 = 1.5f - fx.array();
    Array2f arr_v1 = fx.array() - 1.0f;
    Array2f arr_v2 = fx.array() - 0.5f;
    Array2f ww[3] = {0.5f*arr_v0*arr_v0, 0.75f-arr_v1*arr_v1, 0.5f*arr_v2*arr_v2};

    p.velocity.setZero();
    p.Bp.setZero();

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            Vector2f dpos = Vector2f(i, j) - fx;
            float weight = ww[i][0]*ww[j][1];

            int i2 = i+base_coord_i[0]-gridX_offset;
            int j2 = j+base_coord_i[1];
            int idx_gridnode = j2 + i2*gridY;

            Vector2f node_velocity;
            node_velocity[0] = buffer_grid[1*pitch_grid + idx_gridnode];
            node_velocity[1] = buffer_grid[2*pitch_grid + idx_gridnode];
            p.velocity += weight * node_velocity;
            p.Bp += (4.f*h_inv)*weight *(node_velocity*dpos.transpose());
        }

    // Advection and update of the deformation gradient
    p.pos += p.velocity * dt;
    p.Fe = (Matrix2f::Identity() + dt*p.Bp) * p.Fe;     // p.Bp is the gradient of the velocity vector (it seems)

    ComputePQ(p, kappa, mu);    // pre-computes USV, p, q, etc.
    if(!(p.utility_data & status_crushed)) CheckIfPointIsInsideFailureSurface(p);
    if(p.utility_data & status_crushed)
    {
        ComputeSVD(p, kappa, mu);    // pre-computes USV, p, q, etc.
        Wolper_Drucker_Prager(p);
    }


    // distribute the values of p back into GPU memory: pos, velocity, BP, Fe, Jp_inv, PQ
    for(int i=0; i<icy::SimParams::dim; i++)
    {
        buffer_pts[pt_idx + pitch_pts*(icy::SimParams::posx+i)] = p.pos[i];
        buffer_pts[pt_idx + pitch_pts*(icy::SimParams::velx+i)] = p.velocity[i];
        for(int j=0; j<icy::SimParams::dim; j++)
        {
            buffer_pts[pt_idx + pitch_pts*(icy::SimParams::Fe00 + i*icy::SimParams::dim + j)] = p.Fe(i,j);
            buffer_pts[pt_idx + pitch_pts*(icy::SimParams::Bp00 + i*icy::SimParams::dim + j)] = p.Bp(i,j);
        }
    }

    buffer_pts[pt_idx + pitch_pts*icy::SimParams::idx_Jp_inv] = p.Jp_inv;
    *ptr = p.utility_data; // includes crushed/disable status and grain number

    // at the end of each cycle, PQ are recorded for visualization
    if(recordPQ)
    {
        buffer_pts[pt_idx + pitch_pts*icy::SimParams::idx_P] = p.p_tr;
        buffer_pts[pt_idx + pitch_pts*icy::SimParams::idx_Q] = p.q_tr;
    }
}



// =========================================  DEVICE FUNCTIONS




__device__ void svd(const float a[4], float u[4], float sigma[2], float v[4])
{
    GivensRotation<float> gv(0, 1);
    GivensRotation<float> gu(0, 1);
    singular_value_decomposition(a, gu, sigma, gv);
    gu.template fill<2, float>(u);
    gv.template fill<2, float>(v);
}

__device__ void svd2x2(const Matrix2f &mA, Matrix2f &mU, Vector2f &mS, Matrix2f &mV)
{
    float U[4], V[4], S[2];
    float a[4] = {mA(0,0), mA(0,1), mA(1,0), mA(1,1)};
    svd(a, U, S, V);
    mU << U[0],U[1],U[2],U[3];
    mS << S[0],S[1];
    mV << V[0],V[1],V[2],V[3];
}

__device__ void ComputeSVD(icy::Point &p, const float &kappa, const float &mu)
{
    svd2x2(p.Fe, p.U, p.vSigma, p.V);
    p.vSigmaSquared = p.vSigma.array().square().matrix();
    p.v_s_hat_tr = mu/p.Je_tr * dev_d(p.vSigmaSquared); //mu * pow(Je_tr,-2./d)* dev(SigmaSquared);
}

__device__ void ComputePQ(icy::Point &p, const float &kappa, const float &mu)
{
    Matrix2f &F = p.Fe;
    p.Je_tr = F.determinant();
    p.p_tr = -(kappa/2.f) * (p.Je_tr*p.Je_tr - 1.f);

    p.q_tr = coeff1*mu*(1.f/p.Je_tr)*dev(F*F.transpose()).norm();
//    p.v_s_hat_tr = mu/p.Je_tr * dev_d(p.vSigmaSquared); //mu * pow(Je_tr,-2./d)* dev(SigmaSquared);
}


__device__ void GetParametersForGrain(short grain, float &pmin, float &pmax, float &qmax,
                                      float &beta, float &mSq, float &pmin2)
{
    float var2 = 1.0f + gprms.GrainVariability*0.033f*(-15 + (grain+3)%30);
    float var3 = 1.0f + gprms.GrainVariability*0.1f*(-10 + (grain+4)%11);

    pmax = gprms.IceCompressiveStrength;// * var1;
    pmin = -gprms.IceTensileStrength;// * var2;
    qmax = gprms.IceShearStrength * var3;
    pmin2 = -gprms.IceTensileStrength2 * var2;

    beta = gprms.NACC_beta;
    mSq = (4.f*qmax*qmax*(1.f+2.f*beta))/((pmax-pmin)*(pmax-pmin));
}








// deviatoric part of a diagonal matrix
__device__ Vector2f dev_d(Vector2f Adiag)
{
    return Adiag - Adiag.sum()/2.f*Vector2f::Constant(1.f);
}

__device__ Eigen::Matrix2f dev(Eigen::Matrix2f A)
{
    return A - A.trace()/2*Eigen::Matrix2f::Identity();
}



__device__ void CheckIfPointIsInsideFailureSurface(icy::Point &p)
{
    float beta, M_sq, pmin, pmax, qmax, pmin2;
    GetParametersForGrain(p.grain, pmin, pmax, qmax, beta, M_sq, pmin2);

    if(p.p_tr<0)
    {
        if(p.p_tr<pmin2) {p.utility_data |= status_crushed; return;}
        float q0 = 2*sqrt(-pmax*pmin)*qmax/(pmax-pmin);
        float k = -q0/pmin2;
        float q_limit = k*(p.p_tr-pmin2);
        if(p.q_tr > q_limit) {p.utility_data |= status_crushed; return;}
    }
    else
    {
        float y = (1.f+2.f*beta)*p.q_tr*p.q_tr + M_sq*(p.p_tr + beta*pmax)*(p.p_tr - pmax);
        if(y > 0)
        {
            p.utility_data |= status_crushed;
        }
    }
}



__device__ void Wolper_Drucker_Prager(icy::Point &p)
{
    const float &mu = gprms.mu;
    const float &kappa = gprms.kappa;
    const float &tan_phi = gprms.DP_tan_phi;
    const float &DP_threshold_p = gprms.DP_threshold_p;

    const float &pmax = gprms.IceCompressiveStrength;
    const float &qmax = gprms.IceShearStrength;


    if(p.p_tr < -DP_threshold_p || p.Jp_inv < 1.f)
    {
        // tear in tension or compress until original state
        float p_new = -DP_threshold_p;
        float Je_new = sqrt(-2.f*p_new/kappa + 1.f);
        float sqrt_Je_new = sqrt(Je_new);

        Vector2f vSigma_new(sqrt_Je_new,sqrt_Je_new); //= Vector2d::Constant(1.)*sqrt(Je_new);  //Matrix2d::Identity() * pow(Je_new, 1./(double)d);
        p.Fe = p.U*vSigma_new.asDiagonal()*p.V.transpose();
        p.Jp_inv *= Je_new/p.Je_tr;
    }
    else
    {
        float q_n_1;

        if(p.p_tr > pmax)
        {
            q_n_1 = 0;
        }
        else
        {
            float q_from_dp = (p.p_tr+DP_threshold_p)*tan_phi;
            //q_n_1 = min(q_from_dp,qmax);

            const float pmin = -gprms.IceTensileStrength;
            float q_from_failure_surface = 2.f*sqrt((pmax-p.p_tr)*(p.p_tr-pmin))*qmax/(pmax-pmin);
            q_n_1 = min(q_from_failure_surface, q_from_dp);
        }

        if(p.q_tr >= q_n_1)
        {
            // project onto YS
            float s_hat_n_1_norm = q_n_1*coeff1_inv;
            //Matrix2d B_hat_E_new = s_hat_n_1_norm*(pow(Je_tr,2./d)/mu)*s_hat_tr.normalized() + Matrix2d::Identity()*(SigmaSquared.trace()/d);
            Vector2f vB_hat_E_new = s_hat_n_1_norm*(p.Je_tr/mu)*p.v_s_hat_tr.normalized() +
                                    Vector2f::Constant(1.f)*(p.vSigmaSquared.sum()/d);
            Vector2f vSigma_new = vB_hat_E_new.array().sqrt().matrix();
            p.Fe = p.U*vSigma_new.asDiagonal()*p.V.transpose();
        }
    }
}
