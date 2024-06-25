#include "amr-wind/physics/multiphase/SloshingTank.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/core/field_ops.H"

namespace amr_wind {

SloshingTank::SloshingTank(CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
    , m_pressure(sim.repo().get_field("p"))
{
    amrex::ParmParse pp(identifier());
    pp.query("amplitude", m_amplitude);
    pp.query("peak_enhance", m_kappa);
    pp.query("water_level", m_waterlevel);

    const auto& mphase = sim.physics_manager().get<MultiPhase>();
    m_rho1 = mphase.rho1();
    m_rho2 = mphase.rho2();
    m_gravity = mphase.gravity();
}

/** Initialize the velocity and levelset fields at the beginning of the
 *  simulation. Initializing pressure has little effect for pure AMR-Wind
 *  simulations, but it improves initialization in the overset solver.
 */
void SloshingTank::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    velocity.setVal(0.0);

    auto& levelset = m_levelset(level);
    auto& pressure = m_pressure(level);
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    const auto& probhi = geom.ProbHiArray();
    const amrex::Real Amp = m_amplitude;
    const amrex::Real kappa = m_kappa;
    const amrex::Real water_level = m_waterlevel;
    const amrex::Real Lx = probhi[0] - problo[0];
    const amrex::Real Ly = probhi[1] - problo[1];
    const amrex::Real rho1 = m_rho1;
    const amrex::Real rho2 = m_rho2;
    const amrex::Real grav_z = m_gravity[2];

    for (amrex::MFIter mfi(levelset); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        auto phi = levelset.array(mfi);
        auto p = pressure.array(mfi);

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                const amrex::Real z0 =
                    water_level +
                    Amp * std::exp(
                              -kappa * (std::pow(x - problo[0] - 0.5 * Lx, 2) +
                                        std::pow(y - problo[1] - 0.5 * Ly, 2)));
                phi(i, j, k) = z0 - z;
            });
        amrex::Box const& nbx = mfi.grownnodaltilebox();
        amrex::ParallelFor(
            nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // For pressure nodes, no offset
                const amrex::Real x = problo[0] + i * dx[0];
                const amrex::Real y = problo[1] + j * dx[1];
                const amrex::Real z = problo[2] + k * dx[2];
                const amrex::Real z0 =
                    water_level +
                    Amp * std::exp(
                              -kappa * (std::pow(x - problo[0] - 0.5 * Lx, 2) +
                                        std::pow(y - problo[1] - 0.5 * Ly, 2)));
                // Integrated (top-down in z) phase heights to pressure node
                amrex::Real ih_g =
                    amrex::max(0.0, amrex::min(probhi[2] - z0, probhi[2] - z));
                amrex::Real ih_l =
                    amrex::max(0.0, amrex::min(z0 - z, z0 - problo[2]));
                // Integrated rho at pressure node
                const amrex::Real irho = rho1 * ih_l + rho2 * ih_g;

                // Add term to reference pressure
                p(i, j, k) = -irho * grav_z;
            });
    }
}

} // namespace amr_wind
