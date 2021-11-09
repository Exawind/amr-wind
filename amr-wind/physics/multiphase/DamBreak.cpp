#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/physics/multiphase/DamBreak.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/core/field_ops.H"

namespace amr_wind {

DamBreak::DamBreak(CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
    , m_density(sim.repo().get_field("density"))
{
    amrex::ParmParse pp(identifier());
    pp.queryarr("location", m_loc, 0, AMREX_SPACEDIM);
    pp.query("width", m_width);
    pp.query("height", m_height);
}

/** Initialize the velocity and levelset fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::DamBreakFieldInit
 */
void DamBreak::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    velocity.setVal(0.0, 0, AMREX_SPACEDIM);
    auto& levelset = m_levelset(level);
    auto& density = m_density(level);
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    const amrex::Real xc = m_loc[0];
    const amrex::Real zc = m_loc[2];
    const amrex::Real width = m_width;
    const amrex::Real height = m_height;

    const auto& mphase = m_sim.physics_manager().get<MultiPhase>();
    const amrex::Real rho1 = mphase.rho1();
    const amrex::Real rho2 = mphase.rho2();

    for (amrex::MFIter mfi(levelset); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        auto phi = levelset.array(mfi);
        auto rho = density.array(mfi);
        const amrex::Real eps = std::cbrt(2. * dx[0] * dx[1] * dx[2]);

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                if (x - xc < width && z - zc < height) {
                    phi(i, j, k) = amrex::min(width - x, height - z);
                } else if (x - xc < width && z - zc > height) {
                    phi(i, j, k) = height - (z - zc);
                } else if (x - xc > width && z - zc < height) {
                    phi(i, j, k) = width - (x - xc);
                } else {
                    phi(i, j, k) =
                        amrex::min(width - (x - xc), height - (z - zc));
                }
                amrex::Real smooth_heaviside;
                if (phi(i, j, k) > eps) {
                    smooth_heaviside = 1.0;
                } else if (phi(i, j, k) < -eps) {
                    smooth_heaviside = 0.;
                } else {
                    smooth_heaviside =
                        0.5 *
                        (1.0 + phi(i, j, k) / eps +
                         1.0 / M_PI * std::sin(phi(i, j, k) * M_PI / eps));
                }
                rho(i, j, k) =
                    rho1 * smooth_heaviside + rho2 * (1.0 - smooth_heaviside);
            });
    }
}

} // namespace amr_wind
