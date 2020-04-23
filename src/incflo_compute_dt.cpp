#include <incflo.H>

#include <cmath>
#include <limits>

using namespace amrex;

//
// Compute new dt by using the formula derived in
// "A Boundary Condition Capturing Method for Multiphase Incompressible Flow"
// by Kang et al. (JCP).
//
//  dt/2 * ( C+V + sqrt( (C+V)**2 + 4Fx/dx + 4Fy/dy + 4Fz/dz )
//
// where
//
// C = max(|U|)/dx + max(|V|)/dy + max(|W|)/dz    --> Convection
//
// V = 2 * max(eta/rho) * (1/dx^2 + 1/dy^2 +1/dz^2) --> Diffusion
//
// Fx, Fy, Fz = net acceleration due to external forces
//
// WARNING: We use a slightly modified version of C in the implementation below
//
void incflo::ComputeDt (bool explicit_diffusion)
{
    BL_PROFILE("amr-wind::incflo::ComputeDt")

    Real conv_cfl = 0.0;
    Real diff_cfl = 0.0;
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        auto const dxinv = geom[lev].InvCellSizeArray();
        MultiFab const& vel = velocity()(lev);
        MultiFab const& rho = density()(lev);
        Real conv_lev = 0.0;
        Real diff_lev = 0.0;

        conv_lev = amrex::ReduceMax(vel, 0,
                   [=] AMREX_GPU_HOST_DEVICE (Box const& b,
                                              Array4<Real const> const& v) -> Real
                   {
                       Real mx = -1.0;
                       amrex::Loop(b, [=,&mx] (int i, int j, int k) noexcept
                       {
                           mx = amrex::max(std::abs(v(i,j,k,0))*dxinv[0],
                                           std::abs(v(i,j,k,1))*dxinv[1],
                                           std::abs(v(i,j,k,2))*dxinv[2], mx);
                       });
                       return mx;
                   });
        if (explicit_diffusion) {
            diff_lev = amrex::ReduceMax(rho, 0,
                       [=] AMREX_GPU_HOST_DEVICE (Box const& b,
                                                  Array4<Real const> const& r) -> Real
                       {
                           Real mx = -1.0;
                           amrex::Loop(b, [=,&mx] (int i, int j, int k) noexcept
                           {
                               mx = amrex::max(1.0/r(i,j,k), mx);
                           });
                           return mx;
                       });
            diff_lev *= m_mu;
        }
        
        conv_cfl = std::max(conv_cfl, conv_lev);
        diff_cfl = std::max(diff_cfl, diff_lev*2.*(dxinv[0]*dxinv[0]+dxinv[1]*dxinv[1]+
                                                   dxinv[2]*dxinv[2]));
    }

    Real cd_cfl;
    if (explicit_diffusion) {
        ParallelAllReduce::Max<Real>({conv_cfl,diff_cfl},
                                     ParallelContext::CommunicatorSub());
        cd_cfl = conv_cfl + diff_cfl;
    } else {
        ParallelAllReduce::Max<Real>(conv_cfl,
                                     ParallelContext::CommunicatorSub());
        cd_cfl = conv_cfl;
    }

    // Forcing term
    const auto dxinv_finest = Geom(finest_level).InvCellSizeArray();
    // fixme should we use forces in our dt estimate?
    // abl_godunov_explicit regression test will fail if this is deleted
    Real forc_cfl = std::abs(m_gravity[0] - std::abs(m_gp0[0])) * dxinv_finest[0]
                  + std::abs(m_gravity[1] - std::abs(m_gp0[1])) * dxinv_finest[1]
                  + std::abs(m_gravity[2] - std::abs(m_gp0[2])) * dxinv_finest[2];

    // Combined CFL conditioner
    Real comb_cfl = cd_cfl + std::sqrt(cd_cfl*cd_cfl + 4.0 * forc_cfl);

    m_time.set_current_cfl(comb_cfl);
}
