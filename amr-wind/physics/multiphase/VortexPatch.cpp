#include "amr-wind/physics/multiphase/VortexPatch.H"
#include "amr-wind/physics/multiphase/VortexPatchFieldInit.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

VortexPatch::VortexPatch(CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
{
    // This shouldn't be here, but this is part of the prescirbed velocity field
    // and doesn't fit within VortexPatchFieldInit either.
    amrex::ParmParse pp_vortex_patch("VortexPatch");
    pp_vortex_patch.query("period", m_TT);

    // Instantiate the VortexPatch field initializer
    m_field_init.reset(new VortexPatchFieldInit());
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

    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        (*m_field_init)(vbx, geom, velocity.array(mfi), levelset.array(mfi));
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
            const auto& vbx = mfi.validbox();
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
}

} // namespace amr_wind
