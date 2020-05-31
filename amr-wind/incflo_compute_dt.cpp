#include "amr-wind/incflo.H"

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
void incflo::ComputeDt(bool explicit_diffusion)
{
    BL_PROFILE("amr-wind::incflo::ComputeDt")

    Real conv_cfl = 0.0;
    Real diff_cfl = 0.0;
    Real force_cfl = 0.0;

    const auto& den = density();

    for (int lev = 0; lev <= finest_level; ++lev) {
        auto const dxinv = geom[lev].InvCellSizeArray();
        MultiFab const& vel = icns().fields().field(lev);
        MultiFab const& vel_force = icns().fields().src_term(lev);
        MultiFab const& mu = icns().fields().mueff(lev);
        MultiFab const& rho = den(lev);

        Real conv_lev = 0.0;
        Real diff_lev = 0.0;
        Real force_lev = 0.0;

        conv_lev = amrex::ReduceMax(
            vel, 0,
            [=] AMREX_GPU_HOST_DEVICE(
                Box const& b, Array4<Real const> const& v) -> Real {
                Real mx = -1.0;
                amrex::Loop(b, [=, &mx](int i, int j, int k) noexcept {
                    mx = amrex::max(
                        amrex::Math::abs(v(i, j, k, 0)) * dxinv[0],
                        amrex::Math::abs(v(i, j, k, 1)) * dxinv[1],
                        amrex::Math::abs(v(i, j, k, 2)) * dxinv[2], mx);
                });
                return mx;
            });

        if (explicit_diffusion) {

            const Real dxinv2 =
                2.0 * (dxinv[0] * dxinv[0] + dxinv[1] * dxinv[1] +
                       dxinv[2] * dxinv[2]);

            diff_lev = amrex::ReduceMax(
                rho, mu, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    Box const& b, Array4<Real const> const& rho_arr,
                    Array4<Real const> const& mu_arr) -> Real {
                    Real mx = -1.0;
                    amrex::Loop(b, [=, &mx](int i, int j, int k) noexcept {
                        mx = amrex::max(
                            mu_arr(i, j, k) * dxinv2 / rho_arr(i, j, k), mx);
                    });
                    return mx;
                });
        }

        if (m_time.use_force_cfl()) {
            force_lev = amrex::ReduceMax(
                vel_force, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    Box const& b, Array4<Real const> const& vf) -> Real {
                    Real mx = -1.0;
                    amrex::Loop(b, [=, &mx](int i, int j, int k) noexcept {
                        mx = amrex::max(
                            amrex::Math::abs(vf(i, j, k, 0)) * dxinv[0],
                            amrex::Math::abs(vf(i, j, k, 1)) * dxinv[1],
                            amrex::Math::abs(vf(i, j, k, 2)) * dxinv[2], mx);
                    });
                    return mx;
                });
        }

        conv_cfl = amrex::max(conv_cfl, conv_lev);
        diff_cfl = amrex::max(diff_cfl, diff_lev);
        force_cfl = amrex::max(force_cfl, force_lev);
    }

    ParallelAllReduce::Max<Real>(conv_cfl, ParallelContext::CommunicatorSub());
    if (explicit_diffusion) {
        ParallelAllReduce::Max<Real>(
            diff_cfl, ParallelContext::CommunicatorSub());
    }
    if (m_time.use_force_cfl()) {
        ParallelAllReduce::Max<Real>(
            force_cfl, ParallelContext::CommunicatorSub());
    }

    const Real cd_cfl = conv_cfl + diff_cfl;

    // Combined CFL conditioner
    const Real comb_cfl =
        2.0 * cd_cfl + std::sqrt(cd_cfl * cd_cfl + 4.0 * force_cfl);

    if (m_verbose > 2) {
        amrex::Print() << "conv_cfl: " << conv_cfl << " diff_cfl: " << diff_cfl
                       << " force_cfl: " << force_cfl
                       << " comb_cfl: " << comb_cfl << std::endl;
    }

    m_time.set_current_cfl(comb_cfl);
}
