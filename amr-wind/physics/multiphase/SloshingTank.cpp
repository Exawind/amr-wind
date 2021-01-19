#include "amr-wind/physics/multiphase/SloshingTank.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/core/field_ops.H"

namespace amr_wind {

SloshingTank::SloshingTank(CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
{
    amrex::ParmParse pp(identifier());
    pp.query("amplitude", m_amplitude);
    pp.query("wavelength", m_wavelength);
    pp.query("water_level", m_waterlevel);
}

/** Initialize the velocity and levelset fields at the beginning of the
 *  simulation.
 *
 */
void SloshingTank::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    velocity.setVal(0.0, 0, AMREX_SPACEDIM);

    auto& levelset = m_levelset(level);
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    const amrex::Real Amp = m_amplitude;
    const amrex::Real waveL = m_wavelength;
    const amrex::Real water_level = m_waterlevel;

    for (amrex::MFIter mfi(levelset); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        auto phi = levelset.array(mfi);

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                const amrex::Real z0 =
                    water_level + Amp * std::sin(2.0 * M_PI * x / waveL);
                phi(i, j, k) = z0 - z;
            });
    }
}

} // namespace amr_wind
