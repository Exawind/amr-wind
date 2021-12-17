#include "amr-wind/physics/multiphase/RainDrop.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

RainDrop::RainDrop(CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
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

    auto& mphase = m_sim.physics_manager().get<MultiPhase>();
    // Get VOF array if it is valid
    bool yes_vof = m_sim.repo().field_exists("vof");
    amrex::MultiFab* work_ptr;
    if (yes_vof) {
        auto& volfrac = m_sim.repo().get_field("vof")(level);
        mphase.levelset2vof();
        work_ptr = &volfrac;
    } else {
        work_ptr = &levelset;
    }

    // Use density to calculate mixed velocity in multiphase cells
    amrex::Real rhol = mphase.rho1();
    amrex::Real rhog = mphase.rho2();
    const amrex::Real eps = std::cbrt(2. * dx[0] * dx[1] * dx[2]);
    amrex::Gpu::DeviceVector<amrex::Real> dvel(AMREX_SPACEDIM);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_vel.begin(), m_vel.end(), dvel.begin());
    const auto* vptr = dvel.data();
    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        auto vel = velocity.array(mfi);
        auto work_arr = work_ptr->array(mfi);
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Calculate mass-weighted velocity (gas vel is 0)
                amrex::Real vof;
                if (yes_vof) {
                    // Work array here is vof
                    vof = work_arr(i, j, k);
                } else {
                    // Work array here is phi
                    if (work_arr(i, j, k) > eps) {
                        vof = 1.0;
                    } else if (work_arr(i, j, k) < -eps) {
                        vof = 0.;
                    } else {
                        vof = 0.5 *
                              (1.0 + work_arr(i, j, k) / eps +
                               1.0 / M_PI *
                                   std::sin(work_arr(i, j, k) * M_PI / eps));
                    }
                }
                amrex::Real dens = (1.0 - vof) * rhog + vof * rhol;
                vel(i, j, k, 0) = vof * rhol * vptr[0] / dens;
                vel(i, j, k, 1) = vof * rhol * vptr[1] / dens;
                vel(i, j, k, 2) = vof * rhol * vptr[2] / dens;
            });
    }
}

} // namespace amr_wind
