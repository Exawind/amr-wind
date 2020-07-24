#include "amr-wind/incflo.H"

#include <cmath>
#include <limits>

using namespace amrex;

/** Estimate the new timestep for adaptive timestepping algorithm
 *
 *  \param explicit_diffusion Flag indicating whether user has selected explicit
 *  treatment of diffusion term.
 *
 * Compute new \f$\Delta t\f$ by using the formula derived in "A Boundary Condition
 * Capturing Method for Multiphase Incompressible Flow" by Kang et al. (JCP).
 *
 *  \f{align}
 *  \text{CFL} &= \frac{\Delta t}{2}  \left[\left(C+V \right) + \sqrt{\left(C+V\right)^2 + 4
 * \left(\frac{|F_x|}{\Delta x} + \frac{|F_y|}{\Delta y} + \frac{|F_z|}{\Delta
 * z} \right)} \right] && \\
 *  C &= \text{max} \left(\frac{|U_x|}{\Delta x} ,\ \frac{|U_y|}{\Delta y} ,\ \frac{|U_z|}{\Delta z}\right) && \text{Convection} \\
 *  V &= 2 \left(\frac{\mu}{\rho}\right)_\mathrm{max} \left[\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2} + \frac{1}{\Delta z^2} \right] && \text{Diffusion}
 *  \f}
 *
 *  \f$F_x, F_y, F_z\f$ are acceleration due to forcing terms. By default,
 *  \f$C\f$ is always used for computing CFL. The term \f$V\f$ is only used when
 *  the user has chosen implicit or Crank-Nicholson scheme for diffusion, and
 *  contributions from forcing term when `time.use_force_cfl` is `true` (default
 *  is `true`).
 *
 */
void incflo::ComputeDt(bool explicit_diffusion)
{
    BL_PROFILE("amr-wind::incflo::ComputeDt");

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
