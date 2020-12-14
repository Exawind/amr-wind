#include "amr-wind/physics/multiphase/VortexPatch.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

VortexPatch::VortexPatch(CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
{
    amrex::ParmParse pp(identifier());
    pp.queryarr("location", m_loc, 0, AMREX_SPACEDIM);
    pp.query("radius", m_radius);
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
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();

    const amrex::Real xc = m_loc[0];
    const amrex::Real yc = m_loc[1];
    const amrex::Real zc = m_loc[2];
    const amrex::Real radius = m_radius;
    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        auto vel = velocity.array(mfi);
        auto phi = levelset.array(mfi);
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
