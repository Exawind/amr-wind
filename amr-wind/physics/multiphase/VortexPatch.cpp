#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/physics/multiphase/VortexPatch.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

VortexPatch::VortexPatch(CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
    , m_density(sim.repo().get_field("density"))
{
    amrex::ParmParse pp(identifier());
    pp.queryarr("location", m_loc, 0, AMREX_SPACEDIM);
    pp.query("radius", m_radius);
    pp.query("period", m_TT);
}

/** Initialize the velocity and levelset fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::VortexPatchFieldInit
 */
void VortexPatch::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& levelset = m_levelset(level);
    auto& density = m_density(level);
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();

    const amrex::Real xc = m_loc[0];
    const amrex::Real yc = m_loc[1];
    const amrex::Real zc = m_loc[2];
    const amrex::Real radius = m_radius;

    auto& mphase = m_sim.physics_manager().get<MultiPhase>();
    const amrex::Real rho1 = mphase.rho1();
    const amrex::Real rho2 = mphase.rho2();

    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        auto vel = velocity.array(mfi);
        auto phi = levelset.array(mfi);
        auto rho = density.array(mfi);
        const amrex::Real eps = std::cbrt(2. * dx[0] * dx[1] * dx[2]);

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                vel(i, j, k, 0) =
                    2.0 * std::sin(M_PI * x) * std::sin(M_PI * x) *
                    std::sin(2.0 * M_PI * y) * std::sin(2.0 * M_PI * z);
                vel(i, j, k, 1) = -std::sin(M_PI * y) * std::sin(M_PI * y) *
                                  std::sin(2.0 * M_PI * x) *
                                  std::sin(2.0 * M_PI * z);
                vel(i, j, k, 2) = -std::sin(M_PI * z) * std::sin(M_PI * z) *
                                  std::sin(2.0 * M_PI * x) *
                                  std::sin(2.0 * M_PI * y);

                phi(i, j, k) =
                    radius - std::sqrt(
                                 (x - xc) * (x - xc) + (y - yc) * (y - yc) +
                                 (z - zc) * (z - zc));
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

void VortexPatch::pre_advance_work() {}

void VortexPatch::post_advance_work()
{
    const auto& time = m_sim.time().current_time();

    const int nlevels = m_sim.repo().num_active_levels();
    const auto& geom = m_sim.mesh().Geom();

    // Overriding the velocity field
    for (int lev = 0; lev < nlevels; ++lev) {
        for (amrex::MFIter mfi(m_velocity(lev)); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.growntilebox();
            const auto& dx = geom[lev].CellSizeArray();
            const auto& problo = geom[lev].ProbLoArray();
            const amrex::Real TT = m_TT;
            auto vel = m_velocity(lev).array(mfi);
            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                    const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                    vel(i, j, k, 0) =
                        2.0 * std::sin(M_PI * x) * std::sin(M_PI * x) *
                        std::sin(2.0 * M_PI * y) * std::sin(2.0 * M_PI * z) *
                        std::cos(M_PI * time / TT);
                    vel(i, j, k, 1) = -std::sin(M_PI * y) * std::sin(M_PI * y) *
                                      std::sin(2.0 * M_PI * x) *
                                      std::sin(2.0 * M_PI * z) *
                                      std::cos(M_PI * time / TT);
                    vel(i, j, k, 2) = -std::sin(M_PI * z) * std::sin(M_PI * z) *
                                      std::sin(2.0 * M_PI * x) *
                                      std::sin(2.0 * M_PI * y) *
                                      std::cos(M_PI * time / TT);
                });
        }
    }

    m_velocity.fillpatch(time);
}

} // namespace amr_wind
