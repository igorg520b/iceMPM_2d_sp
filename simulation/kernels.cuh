#include "parameters_sim.h"


using namespace Eigen;

constexpr t_PointReal d = 2; // dimensions
constexpr t_PointReal coeff1 = 1.4142135623730950; // sqrt((6-d)/2.);
constexpr t_PointReal coeff1_inv = 0.7071067811865475;
constexpr uint32_t status_crushed = 0x10000;
constexpr uint32_t status_disabled = 0x20000;

__device__ uint8_t gpu_error_indicator;
__device__ int gpu_disabled_points_count;
__constant__ SimParams gprms;

// device-side parameters that are typically reused in kernels
__device__ t_PointReal *dev_buffer_pts;


__global__ void partition_kernel_p2g(const int gridX, const int pitch_grid,
                                     const int count_pts, const int pitch_pts,
                                     t_GridReal *buffer_grid)
{
    int pt_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(pt_idx >= count_pts) return;

    uint32_t utility_data = *reinterpret_cast<const uint32_t*>(&dev_buffer_pts[pt_idx + pitch_pts*SimParams::idx_utility_data]);
    if(utility_data & status_disabled) return; // point is disabled

    const double &h = gprms.cellsize;
    const int &gridY = gprms.GridYTotal;

    // pull point data from SOA
    PointVector2r pos, velocity;
    PointMatrix2r Bp, Fe;

    for(int i=0; i<SimParams::dim; i++)
    {
        pos[i] = dev_buffer_pts[pt_idx + pitch_pts*(SimParams::posx+i)];
        velocity[i] = dev_buffer_pts[pt_idx + pitch_pts*(SimParams::velx+i)];
        for(int j=0; j<SimParams::dim; j++)
        {
            Fe(i,j) = dev_buffer_pts[pt_idx + pitch_pts*(SimParams::Fe00 + i*SimParams::dim + j)];
            Bp(i,j) = dev_buffer_pts[pt_idx + pitch_pts*(SimParams::Bp00 + i*SimParams::dim + j)];
        }
    }
    t_PointReal Jp_inv = dev_buffer_pts[pt_idx + pitch_pts*SimParams::idx_Jp_inv];
    const t_PointReal thickness = dev_buffer_pts[pt_idx + pitch_pts*SimParams::idx_thickness];
    double particle_mass = gprms.ParticleMass * thickness;

    const uint32_t cell = *reinterpret_cast<const uint32_t*>(&dev_buffer_pts[pt_idx + pitch_pts*SimParams::integer_cell_idx]);
    Eigen::Vector2i cell_i((int)(cell & 0xffff), (int)(cell >> 16));

    PointMatrix2r PFt;
    PFt = KirchhoffStress_Wolper(Fe);
    // PFt = Water(buffer_pts[pt_idx + pitch_pts*SimParams::idx_Jp_inv]);
//    PointMatrix2r subterm2 = particle_mass*Bp - (gprms.dt_vol_Dpinv)*PFt;

    // version that accounts for surface density change
    PointMatrix2r subterm2 = particle_mass*Bp - (gprms.dt_vol_Dpinv*Jp_inv*thickness)*PFt;

    PointArray2r ww[3];
    CalculateWeightCoeffs(pos, ww);

    for (int i = -1; i <= 1; i++)
        for (int j = -1; j <= 1; j++)
        {
            t_PointReal Wip = ww[i+1][0]*ww[j+1][1];
            t_PointReal incM = Wip*particle_mass;
            PointVector2r dpos((i-pos[0])*h, (j-pos[1])*h);
            PointVector2r incV = Wip*(velocity*particle_mass + subterm2*dpos);

            // the x-index of the cell takes into accout the partition's offset of the gird fragment
            int idx_gridnode = (j+cell_i[1]) + (i+cell_i[0])*gridY;

            // Udpate mass, velocity and force
            atomicAdd(&buffer_grid[0*pitch_grid + idx_gridnode], (t_GridReal)incM);
            atomicAdd(&buffer_grid[1*pitch_grid + idx_gridnode], (t_GridReal)incV[0]);
            atomicAdd(&buffer_grid[2*pitch_grid + idx_gridnode], (t_GridReal)incV[1]);
        }

    // check if a point is out of bounds of the grid box
    if(cell_i[0] < 1) gpu_error_indicator = 70;
    if(cell_i[1] < 1) gpu_error_indicator = 71;
    if(cell_i[0] > (gridX-2)) gpu_error_indicator = 72;
    if(cell_i[1] > gridY-2) gpu_error_indicator = 73;
}




__global__ void partition_kernel_update_nodes(const int nNodes, const int pitch_grid,
                                              t_GridReal *buffer_grid,
                                              t_PointReal simulation_time, const uint8_t *grid_status,
                                              GridVector2r vWind, const float interpolation_coeff_w)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= nNodes) return;

    const int &gridY = gprms.GridYTotal;

    t_GridReal mass = buffer_grid[SimParams::grid_idx_mass*pitch_grid + idx];
    if(mass == 0) return;

    t_GridReal vx = buffer_grid[SimParams::grid_idx_px*pitch_grid + idx];
    t_GridReal vy = buffer_grid[SimParams::grid_idx_py*pitch_grid + idx];

    // const double &dt = gprms.InitialTimeStep;
    const double &cellsize = gprms.cellsize;
    const double &vmax = gprms.vmax;

    Vector2i gi(idx/gridY, idx%gridY);   // integer x-y index of the grid node
    GridVector2r gnpos = gi.cast<t_GridReal>()*cellsize;    // position of the grid node in the whole grid

    GridVector2r velocity(vx, vy);
    velocity /= mass;

    uint8_t is_modeled_area = grid_status[idx];
    if(is_modeled_area != 100)
    {
        velocity.setZero();
    }
    else
    {
        //grid_water_current

        t_GridReal vcx = buffer_grid[SimParams::grid_idx_current_vx*pitch_grid + idx];
        t_GridReal vcy = buffer_grid[SimParams::grid_idx_current_vy*pitch_grid + idx];

        GridVector2r wvel(vcx, vcy);
        wvel *= (1+min(simulation_time/(3600*10), 3.));

//        wvel.x() = 0.1;
//        wvel.y() = -0.3;

        // quadratic term
/*
        GridVector2r vrel = wvel - velocity;
        const t_GridReal Density_coeff_A = gprms.currentDragCoeff_waterDensity * hsq;
        t_GridReal dragForce = vrel.squaredNorm() * Density_coeff_A;
        vrel.normalize();
//        vrel *= dragForce * dt / mass;
        vrel *= dragForce * dt;

        const double vmax2 = vmax*0.01;
        if(vrel.norm() > vmax2) vrel = vrel.normalized()*vmax2;

        velocity += vrel;
*/

        // linear term
//        const double coeff = 0.005;
        const double coeff = 0.05;
        velocity = (1-coeff)*velocity + coeff*wvel;

        /*
        GridVector2r wvel2 = wvel-velocity;
        wvel2 = wvel2 * wvel2.norm();
        const double coeff2 = 0.01;
        velocity += coeff2*wvel2;

        GridVector2r wvel3 = (wvel-velocity)*1000*hsq*dt/mass;
        //wvel3 = wvel3 * wvel3.norm();
        if(wvel3.norm() > vmax*0.1) wvel3 = wvel3.normalized()*vmax*0.1;
        const double coeff3 = 0.01;
        velocity += coeff3*wvel3;
*/


        //interpolation_coeff_w

//        float lat = gprms.LatMin + (gprms.LatMax-gprms.LatMin)*(float)gi.y()/(float)gprms.GridY;
//        float lon = gprms.LonMin + (gprms.LonMax-gprms.LonMin)*(float)gi.x()/(float)gprms.GridXTotal;
//        vWind = get_wind_vector(lat, lon, interpolation_coeff);


/*
        // wind drag
        //t_GridReal hsq = gprms.cellsize * gprms.cellsize;
        vWind *= 10 * wvel.norm();
        vrel = vWind - velocity;
        const t_GridReal airDensity_coeff_A = gprms.windDragCoeff_airDensity * hsq;
        dragForce = vrel.squaredNorm() * airDensity_coeff_A;
        vrel.normalize();
        vrel *= dragForce * dt / mass;
        velocity = (1-coeff2)*velocity + coeff2*(velocity+vrel);
*/

    }

    // write the updated grid velocity back to memory
    if(velocity.squaredNorm() > vmax*vmax*0.1*0.1) velocity.setZero();
    buffer_grid[SimParams::grid_idx_px*pitch_grid + idx] = (t_GridReal)velocity[0];
    buffer_grid[SimParams::grid_idx_py*pitch_grid + idx] = (t_GridReal)velocity[1];
}



__global__ void partition_kernel_g2p(const bool recordPQ, const int pitch_grid,
                                     const int count_pts, const int pitch_pts,
                                     const t_GridReal *buffer_grid,
                                     int applyGlensLaw)
{
    const int pt_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(pt_idx >= count_pts) return;

    // skip if a point is disabled
    uint32_t utility_data = *reinterpret_cast<uint32_t*>(&dev_buffer_pts[pt_idx + pitch_pts*SimParams::idx_utility_data]);
    if(utility_data & status_disabled) return; // point is disabled

    const double &h_inv = gprms.cellsize_inv;
    const double &dt = gprms.InitialTimeStep;
    const double &mu = gprms.mu;
    const double &kappa = gprms.kappa;
    const int &gridY = gprms.GridYTotal;
    const int &gridX = gprms.GridXTotal;

    PointVector2r pos;
    PointMatrix2r Fe;

    // pull point data from SOA
    for(int i=0; i<SimParams::dim; i++)
    {
        pos[i] = dev_buffer_pts[pt_idx + pitch_pts*(SimParams::posx+i)];
        for(int j=0; j<SimParams::dim; j++)
        {
            Fe(i,j) = dev_buffer_pts[pt_idx + pitch_pts*(SimParams::Fe00 + i*SimParams::dim + j)];
        }
    }
    t_PointReal initial_thickness = dev_buffer_pts[pt_idx + pitch_pts*SimParams::idx_thickness];
    t_PointReal Jp_inv = dev_buffer_pts[pt_idx + pitch_pts*SimParams::idx_Jp_inv];
    t_PointReal qp = dev_buffer_pts[pt_idx + pitch_pts*SimParams::idx_Qp];
    uint16_t grain = (uint16_t) (utility_data & 0xffff);

    uint32_t cell = *reinterpret_cast<const uint32_t*>(&dev_buffer_pts[pt_idx + pitch_pts*SimParams::integer_cell_idx]);
    // coords of grid node for point
    Eigen::Vector2i cell_i((int)(cell & 0xffff), (int)(cell >> 16));

    // optimized method of computing the quadratic weight function without conditional operators
    PointArray2r ww[3];
    CalculateWeightCoeffs(pos, ww);

    PointVector2r p_velocity;
    PointMatrix2r p_Bp;
    p_velocity.setZero(); p_Bp.setZero();

    for (int i = -1; i <= 1; i++)
        for (int j = -1; j <= 1; j++)
        {
            PointVector2r dpos = PointVector2r(i, j) - pos;
            t_PointReal weight = ww[i+1][0]*ww[j+1][1];

            int idx_gridnode = (j+cell_i[1]) + (i+cell_i[0])*gridY; // grid node index within the 3x3 loop

            PointVector2r node_velocity;
            node_velocity[0] = (t_PointReal)buffer_grid[1*pitch_grid + idx_gridnode];
            node_velocity[1] = (t_PointReal)buffer_grid[2*pitch_grid + idx_gridnode];
            p_velocity += weight * node_velocity;
            p_Bp += (4.f*h_inv*weight) * (node_velocity*dpos.transpose());
        }

    // Advection and update of the deformation gradient
    pos += p_velocity * (dt*h_inv); // position is in local cell coordinates [-0.5 to 0.5]

    // encode the position of the point as coordinates + cell index
    bool cell_updated = false;
    if(pos.x() > 0.5) { pos.x() -= 1.0; cell_i.x()++; cell_updated = true; }
    else if(pos.x() < -0.5) { pos.x() += 1.0; cell_i.x()--; cell_updated = true; }
    if(pos.y() > 0.5) { pos.y() -= 1.0; cell_i.y()++; cell_updated = true; }
    else if(pos.y() < -0.5) { pos.y() += 1.0; cell_i.y()--; cell_updated = true; }

    if(cell_updated)
    {
        if(cell_i.x() <= 1 || cell_i.x() >= gridX-2 || cell_i.y() <= 1 || cell_i.y() >= gridY-2)
        {
            utility_data |= status_disabled;
            atomicAdd(&gpu_disabled_points_count, 1);
        }
    }

    // ensure the coordinates are valid
    if(isnan(pos.x()) || isnan(pos.y())) { gpu_error_indicator = 1; printf("kernel error 1\n");}
    if(pos.x() > 0.5 || pos.x() < -0.5 || pos.y() > 0.5 || pos.y() < -0.5) { gpu_error_indicator = 2; printf("kernel error 2\n"); }

    Fe = (PointMatrix2r::Identity() + dt*p_Bp) * Fe;     // p.Bp is the gradient of the velocity vector (it seems)

    t_PointReal Je_tr, p_tr, q_tr;
    ComputePQ(Je_tr, p_tr, q_tr, kappa, mu, Fe);    // computes P, Q, J

    PointMatrix2r U, V;
    PointVector2r vSigma, vSigmaSquared, v_s_hat_tr;
    ComputeSVD(Fe, U, vSigma, V, vSigmaSquared, v_s_hat_tr, kappa, mu, Je_tr);

    if(!(utility_data & status_crushed)) CheckIfPointIsInsideFailureSurface(utility_data, grain, p_tr, q_tr, initial_thickness);
    if(utility_data & status_crushed)
    {
        Wolper_Drucker_Prager(initial_thickness, p_tr, q_tr, Je_tr, U, V, vSigmaSquared, v_s_hat_tr, Fe, Jp_inv);
    }
//    else if(applyGlensLaw)
//    {
//        t_PointReal dt_special = dt*applyGlensLaw;
//        Glen_Nye_flow_law(dt_special, q_tr, vSigmaSquared, U, V, v_s_hat_tr, Fe, qp);
//    }


    // distribute the values of p back into GPU memory: pos, velocity, BP, Fe, Jp_inv, PQ
    for(int i=0; i<SimParams::dim; i++)
    {
        dev_buffer_pts[pt_idx + pitch_pts*(SimParams::posx+i)] = pos[i];
        dev_buffer_pts[pt_idx + pitch_pts*(SimParams::velx+i)] = p_velocity[i];
        for(int j=0; j<SimParams::dim; j++)
        {
            dev_buffer_pts[pt_idx + pitch_pts*(SimParams::Fe00 + i*SimParams::dim + j)] = Fe(i,j);
            dev_buffer_pts[pt_idx + pitch_pts*(SimParams::Bp00 + i*SimParams::dim + j)] = p_Bp(i,j);
        }
    }

    dev_buffer_pts[pt_idx + pitch_pts*SimParams::idx_Jp_inv] = Jp_inv;
    dev_buffer_pts[pt_idx + pitch_pts*SimParams::idx_Qp] = qp;

    // save crushed/disabled status
    *reinterpret_cast<uint32_t*>(&dev_buffer_pts[pt_idx + pitch_pts*SimParams::idx_utility_data]) = utility_data;

    if(cell_updated)
    {
        cell = ((uint32_t)cell_i[1] << 16) | (uint32_t)cell_i[0];
        *reinterpret_cast<uint32_t*>(&dev_buffer_pts[pt_idx + pitch_pts*SimParams::integer_cell_idx]) = cell;
    }

    // upon request, PQ are recorded for visualization
    if(recordPQ)
    {
        dev_buffer_pts[pt_idx + pitch_pts*SimParams::idx_P] = p_tr;
        dev_buffer_pts[pt_idx + pitch_pts*SimParams::idx_Q] = q_tr;
    }
}



// =========================================  DEVICE FUNCTIONS



__forceinline__ __device__ void CalculateWeightCoeffs(const PointVector2r &pos, PointArray2r ww[3])
{
    // optimized method of computing the quadratic (!) weight function (no conditional operators)
    PointArray2r arr_v0 = 0.5 - pos.array();
    PointArray2r arr_v1 = pos.array();
    PointArray2r arr_v2 = pos.array() + 0.5;
    ww[0] = 0.5*arr_v0*arr_v0;
    ww[1] = 0.75-arr_v1*arr_v1;
    ww[2] = 0.5*arr_v2*arr_v2;
}




__device__ void ComputePQ(t_PointReal &Je_tr, t_PointReal &p_tr, t_PointReal &q_tr,
    const double &kappa, const double &mu, const PointMatrix2r &F)
{
    Je_tr = F.determinant();
    p_tr = -(kappa/2.) * (Je_tr*Je_tr - 1.);
    q_tr = coeff1*mu*(1./Je_tr)*dev(F*F.transpose()).norm();
}











__device__ void CheckIfPointIsInsideFailureSurface(uint32_t &utility_data, const uint16_t &grain,
                            const t_PointReal &p, const t_PointReal &q, const t_PointReal &strength)
{
    t_PointReal beta, M_sq, pmin, pmax, qmax, pmin2;
    GetParametersForGrain(grain, pmin, pmax, qmax, beta, M_sq, pmin2);
//    qmax *= strength;

    if(p<0)
    {
        if(p<pmin2) { utility_data |= status_crushed; return; }
        t_PointReal q0 = 2*sqrt(-pmax*pmin)*qmax/(pmax-pmin);
        t_PointReal k = -q0/pmin2;
        t_PointReal q_limit = k*(p-pmin2);
        if(q > q_limit) { utility_data |= status_crushed; return; }
    }
    else
    {
        t_PointReal y = (1.+2.*beta)*q*q + M_sq*(p+beta*pmax) * (p-pmax);
        if(y > 0) utility_data |= status_crushed;
    }
}


__device__ t_PointReal smoothstep(t_PointReal x)
{
    if(x<0) x = 0;
    if(x>1) x = 1;
    return (x*x)*(3.0 - 2.0 * x);
}




__device__ void Wolper_Drucker_Prager(const t_PointReal &initial_thickness,
                                      const t_PointReal &p_tr, const t_PointReal &q_tr, const t_PointReal &Je_tr,
const PointMatrix2r &U, const PointMatrix2r &V, const PointVector2r &vSigmaSquared, const PointVector2r &v_s_hat_tr,
                                      PointMatrix2r &Fe, t_PointReal &Jp_inv)
{
    const double &mu = gprms.mu;
    const double &kappa = gprms.kappa;
    t_PointReal DP_threshold_p = gprms.DP_threshold_p;
//    DP_threshold_p *= Jp_inv;

    const t_PointReal pmin = -gprms.IceTensileStrength;

    const double &pmax = gprms.IceCompressiveStrength;
    const double &qmax = gprms.IceShearStrength;

    const t_PointReal tan_phi = tan(gprms.DP_phi*SimParams::pi/180);

    if(p_tr < DP_threshold_p)
    {
        if(Jp_inv < 0.1)
        {
            t_PointReal sqrt_Je_new = sqrt(Je_tr);
            PointVector2r vSigma_new(sqrt_Je_new,sqrt_Je_new); //= Vector2d::Constant(1.)*sqrt(Je_new);  //Matrix2d::Identity() * pow(Je_new, 1./(double)d);
            Fe = U*vSigma_new.asDiagonal()*V.transpose();
        }
        else
        {
        // stretching in tension - no resistance
            PointVector2r vSigma_new(1.0,1.0);
            Fe = U*vSigma_new.asDiagonal()*V.transpose();
            Jp_inv /= Je_tr;
        }
        return;
    }

    // determine q_yeld from the combination of DP / elliptic yield surface, whichever is lower
    t_PointReal q_yield = 0;
    if(p_tr < pmax)
    {
        t_PointReal q_from_failure_surface = 2*sqrt((pmax-p_tr)*(p_tr-pmin))*qmax/(pmax-pmin);

        t_PointReal q_from_dp = max((double)0, (p_tr-DP_threshold_p)*tan_phi);
        q_yield = min(q_from_failure_surface, q_from_dp);
    }
    else
    {
        // such hight pressures should not happen - everythigng is liquified
    }

    if(q_tr > q_yield)
    {
        // plasticity will be applied

        // estimate the new P based on the ridge height
//        t_PointReal p_ridge_max = SimParams::g * gprms.IceDensity * initial_thickness * Jp_inv;

        t_PointReal p_n_1 = p_tr;





        // re-evaluate q
        const t_PointReal q_from_dp = max((double)0, (p_n_1-DP_threshold_p)*tan_phi);
        const t_PointReal q_from_failure_surface = 2*sqrt((pmax-p_n_1)*(p_n_1-pmin))*qmax/(pmax-pmin);
        const t_PointReal q_n_1 = min(q_from_failure_surface, q_from_dp);

        // given p_n_1 and q_n_1, compute Fe_new
        const t_PointReal Je_new = sqrt(-2*p_n_1/kappa + 1);
        t_PointReal s_hat_n_1_norm = q_n_1*coeff1_inv;
        //Matrix2d B_hat_E_new = s_hat_n_1_norm*(pow(Je_tr,2./d)/mu)*s_hat_tr.normalized() + Matrix2d::Identity()*(SigmaSquared.trace()/d);

        const PointVector2r vB_hat_E_new = s_hat_n_1_norm*(Je_tr/mu)*v_s_hat_tr.normalized() +
                                     PointVector2r::Constant(1)*(Je_new);

        const PointVector2r vSigma_new = vB_hat_E_new.array().sqrt().matrix();
        Fe = U*vSigma_new.asDiagonal()*V.transpose();
        Jp_inv *= Je_new/Je_tr;
    }

//    t_PointReal smstp = smoothstep((Jp_inv-0.5)/2.);


    /*
    if(p_tr < DP_threshold_p || Jp_inv < 1)
    {
        // stretching in tension
        t_PointReal p_new = 0;
        t_PointReal Je_new = sqrt(-2*p_new/kappa + 1);
        t_PointReal sqrt_Je_new = sqrt(Je_new);

        PointVector2r vSigma_new(sqrt_Je_new,sqrt_Je_new); //= Vector2d::Constant(1.)*sqrt(Je_new);  //Matrix2d::Identity() * pow(Je_new, 1./(double)d);
        Fe = U*vSigma_new.asDiagonal()*V.transpose();
        Jp_inv *= Je_new/Je_tr;
    }
    else
    {
        t_PointReal q_n_1;

        if(p_tr > pmax)
        {
            q_n_1 = 0;
        }
        else
        {
            t_PointReal q_from_dp = (p_tr-DP_threshold_p)*tan_phi;
            //q_n_1 = min(q_from_dp,qmax);

            const t_PointReal pmin = -gprms.IceTensileStrength;
            t_PointReal q_from_failure_surface = 2*sqrt((pmax-p_tr)*(p_tr-pmin))*qmax/(pmax-pmin);
            q_n_1 = min(q_from_failure_surface, q_from_dp);
        }

        if(q_tr >= q_n_1)
        {
            // project onto YS
            t_PointReal s_hat_n_1_norm = q_n_1*coeff1_inv;
            //Matrix2d B_hat_E_new = s_hat_n_1_norm*(pow(Je_tr,2./d)/mu)*s_hat_tr.normalized() + Matrix2d::Identity()*(SigmaSquared.trace()/d);
            PointVector2r vB_hat_E_new = s_hat_n_1_norm*(Je_tr/mu)*v_s_hat_tr.normalized() +
                                    PointVector2r::Constant(1)*(vSigmaSquared.sum()/d);
            PointVector2r vSigma_new = vB_hat_E_new.array().sqrt().matrix();
            Fe = U*vSigma_new.asDiagonal()*V.transpose();
        }
    }

*/
}


__device__ void svd2x2(const PointMatrix2r &mA, PointMatrix2r &mU, PointVector2r &mS, PointMatrix2r &mV)
{
    t_PointReal U[4], V[4], S[2];

    GivensRotation<t_PointReal> gv(0, 1);
    GivensRotation<t_PointReal> gu(0, 1);
    singular_value_decomposition(mA.data(), gu, S, gv);
    gu.template fill<2, t_PointReal>(U);
    gv.template fill<2, t_PointReal>(V);

    mU << U[0],U[1],U[2],U[3];
    mS << S[0],S[1];
    mV << V[0],V[1],V[2],V[3];
}

__device__ void ComputeSVD(const PointMatrix2r &Fe, PointMatrix2r &U, PointVector2r &vSigma, PointMatrix2r &V,
                            PointVector2r &vSigmaSquared, PointVector2r &v_s_hat_tr,
                            const t_PointReal &kappa, const t_PointReal &mu, const t_PointReal &Je_tr)
{
    svd2x2(Fe, U, vSigma, V);
    vSigmaSquared = vSigma.array().square().matrix();
    v_s_hat_tr = mu/Je_tr * dev_d(vSigmaSquared); //mu * pow(Je_tr,-2./d)* dev(SigmaSquared);
}



// deviatoric part of a diagonal matrix
__device__ PointVector2r dev_d(PointVector2r Adiag)
{
    return Adiag - Adiag.sum()/2.*PointVector2r::Constant(1.);
}

__device__ PointMatrix2r dev(PointMatrix2r A)
{
    return A - A.trace()/2*PointMatrix2r::Identity();
}



__device__ PointMatrix2r KirchhoffStress_Wolper(const PointMatrix2r &F)
{
    const double &kappa = gprms.kappa;
    const double &mu = gprms.mu;

    // Kirchhoff stress as per Wolper (2019)
    t_PointReal Je = F.determinant();
    PointMatrix2r b = F*F.transpose();
    PointMatrix2r PFt = mu*(1./Je)*dev(b) + kappa*0.5*(Je*Je-1.)*PointMatrix2r::Identity();
    return PFt;
}


__device__ void GetParametersForGrain(uint32_t utility_data, t_PointReal &pmin, t_PointReal &pmax, t_PointReal &qmax,
                                      t_PointReal &beta, t_PointReal &mSq, t_PointReal &pmin2)
{
    //    t_PointReal var2 = 1.0 + gprms.GrainVariability*0.033*(-15 + ((int)grain+3)%30);
    //    t_PointReal var3 = 1.0 + gprms.GrainVariability*0.1*(-10 + ((int)grain+4)%11);

//    bool is_weakened = utility_data & status_weakened;

//    t_PointReal gv = gprms.GrainVariability;

    t_PointReal var2 = 1.0;
    t_PointReal var3 = 1.0;
//    if(((int)grain)%3==0) var2 = 1.0 - gv;
//    if(((int)grain+1)%5==0) var3 = 1.0 - gv;
//    if(is_weakened)
//    {
//        var3 = var2 = 1.0 - gv;
//    }


    pmax = gprms.IceCompressiveStrength;// * var1;
    pmin = -gprms.IceTensileStrength;// * var2;

    qmax = gprms.IceShearStrength * var3;
    pmin2 = -gprms.IceTensileStrength2 * var2;

    beta = gprms.IceTensileStrength/gprms.IceCompressiveStrength;
    mSq = (4.*qmax*qmax*(1.+2.*beta))/((pmax-pmin)*(pmax-pmin));
}

/*
__device__ GridVector2r get_wind_vector(float lat, float lon, float tb)
{
    const double &gridLatMin = gprms.gridLatMin;
    const double &gridLonMin = gprms.gridLonMin;

    // space
    int lat_cell = (int)((lat-gridLatMin)/WindInterpolator::gridCellSize);
    int lon_cell = (int)((lon-gridLonMin)/WindInterpolator::gridCellSize);

    // Compute local coordinates within the cell
    float localLon = lon - (gridLonMin + lon_cell * WindInterpolator::gridCellSize);
    float localLat = lat - (gridLatMin + lat_cell * WindInterpolator::gridCellSize);

    // Compute barycentric coordinates
    float ub = localLon / WindInterpolator::gridCellSize;
    float vb = localLat / WindInterpolator::gridCellSize;

    GridVector2r cell_values0[2][2], cell_values1[2][2];
    for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
        {
            cell_values0[i][j] = GridVector2r(wgrid[lat_cell+i][lon_cell+j][0], wgrid[lat_cell+i][lon_cell+j][1]);
            cell_values1[i][j] = GridVector2r(wgrid[lat_cell+i][lon_cell+j][2], wgrid[lat_cell+i][lon_cell+j][3]);
        }
    GridVector2r ipVal[2];

    ipVal[0] =
        (1 - ub) * (1 - vb) * cell_values0[0][0] +
        ub * (1 - vb) * cell_values0[0][1] +
        (1 - ub) * vb * cell_values0[1][0] +
        ub * vb * cell_values0[1][1];

    ipVal[1] =
        (1 - ub) * (1 - vb) * cell_values1[0][0] +
        ub * (1 - vb) * cell_values1[0][1] +
        (1 - ub) * vb * cell_values1[1][0] +
        ub * vb * cell_values1[1][1];

    GridVector2r final_result = (1-tb)*ipVal[0] + tb*ipVal[1];
    return final_result;
}
*/


/*
__device__ void Glen_Nye_flow_law(const t_PointReal dt, const t_PointReal &q_tr,
const PointVector2r &vSigmaSquared,
const PointMatrix2r &U,
const PointMatrix2r &V,
const PointVector2r &v_s_hat_tr,
                                  PointMatrix2r &Fe, t_PointReal &qp)

{
    const t_PointReal &mu = gprms.mu;
    const t_PointReal &A = gprms.GlenA;


    const t_PointReal Je_tr = Fe.determinant();
    t_PointReal epsilon_dot_dt = A * q_tr*q_tr*q_tr * dt;      // Glen's Law


    t_PointReal q_n_1 = max(q_tr - mu*epsilon_dot_dt, 0.);


    t_PointReal s_hat_n_1_norm = q_n_1*coeff1_inv;
    PointVector2r vB_hat_E_new = s_hat_n_1_norm*(Je_tr/mu)*v_s_hat_tr.normalized() +
                                 PointVector2r::Constant(1)*(vSigmaSquared.sum()/d);
    PointVector2r vSigma_new = vB_hat_E_new.array().sqrt().matrix();
    Fe = U*vSigma_new.asDiagonal()*V.transpose();
    qp *= (q_n_1/q_tr);
}
*/

