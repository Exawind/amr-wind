#include "amr-wind/physics/multiphase/RainDrop.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

RainDrop::RainDrop(CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
    , m_vof(sim.repo().get_field("vof"))
    , m_sim(sim)
{
    amrex::ParmParse pp(identifier());
    pp.queryarr("location", m_loc, 0, AMREX_SPACEDIM);
    pp.query("radius", m_radius);
    pp.queryarr("droplet_velocity", m_vel, 0, AMREX_SPACEDIM);
}

/** Initialize the velocity and levelset fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::RainDropFieldInit
 */
void RainDrop::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& levelset = m_levelset(level);
    auto& volfrac = m_vof(level);
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();

    const amrex::Real xc = m_loc[0];
    const amrex::Real yc = m_loc[1];
    const amrex::Real zc = m_loc[2];
    const amrex::Real radius = m_radius;
    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        auto phi = levelset.array(mfi);
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                // Calculate level-set
                phi(i, j, k) =
                    radius - std::sqrt(
                                 (x - xc) * (x - xc) + (y - yc) * (y - yc) +
                                 (z - zc) * (z - zc));
            });
    }

    // Get VOF array
    auto& mphase = m_sim.physics_manager().get<MultiPhase>();
    mphase.levelset2vof();

    // Use density to calculate mixed velocity in multiphase cells
    amrex::Real rhol = mphase.rho1();
    amrex::Real rhog = mphase.rho2();
    amrex::Gpu::DeviceVector<amrex::Real> dvel(AMREX_SPACEDIM);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_vel.begin(), m_vel.end(), dvel.begin());
    const auto* vptr = dvel.data();
    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        auto vof_arr = volfrac.array(mfi);
        auto vel = velocity.array(mfi);
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Calculate mass-weighted velocity (gas vel is 0)
                amrex::Real vof = vof_arr(i, j, k);
                amrex::Real dens = (1.0 - vof) * rhog + vof * rhol;
                vel(i, j, k, 0) = vof * rhol * vptr[0] / dens;
                vel(i, j, k, 1) = vof * rhol * vptr[1] / dens;
                vel(i, j, k, 2) = vof * rhol * vptr[2] / dens;
            });
    }
}

} // namespace amr_wind
