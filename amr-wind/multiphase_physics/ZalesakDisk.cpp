#include "amr-wind/multiphase_physics/ZalesakDisk.H"
#include "amr-wind/multiphase_physics/ZalesakDiskFieldInit.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/core/field_ops.H"

namespace amr_wind {

ZalesakDisk::ZalesakDisk(CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{
    // This shouldn't be here, but this is part of the prescirbed velocity field
    // and doesn't fit within ZalesakDiskFieldInit either.
    amrex::ParmParse pp_vortex_patch("ZalesakDisk");
    pp_vortex_patch.query("period", m_TT);

    // Register levelset equation
    auto& levelset_eqn = sim.pde_manager().register_transport_pde("Levelset");

    // Defer getting levelset field until PDE has been registered
    m_levelset = &(levelset_eqn.fields().field);

    // Instantiate the ZalesakDisk field initializer
    m_field_init.reset(new ZalesakDiskFieldInit());
}

/** Initialize the velocity and levelset fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::ZalesakDiskFieldInit
 */
void ZalesakDisk::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& density = m_density(level);
    auto& levelset = (*m_levelset)(level);

    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        (*m_field_init)(
            vbx, geom, velocity.array(mfi), density.array(mfi),
            levelset.array(mfi));
    }
}

void ZalesakDisk::pre_advance_work() {}

void ZalesakDisk::post_advance_work()
{

    const int nlevels = m_sim.repo().num_active_levels();
    const auto& geom = m_sim.mesh().Geom();

    // Overriding the velocity field
    for (int lev = 0; lev < nlevels; ++lev) {
        for (amrex::MFIter mfi(m_velocity(lev)); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            const auto& dx = geom[lev].CellSizeArray();
            const auto& problo = geom[lev].ProbLoArray();

            auto vel = m_velocity(lev).array(mfi);
            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];

                    vel(i, j, k, 0) = 2.0 * M_PI / m_TT * (0.5 - y);
                    vel(i, j, k, 1) = 2.0 * M_PI / m_TT * (x - 0.5);
                    vel(i, j, k, 2) = 0.0;
                });
        }
    }
}

} // namespace amr_wind
