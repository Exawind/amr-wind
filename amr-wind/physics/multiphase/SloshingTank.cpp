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
    pp.query("depth", m_depth);
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
    const amrex::Real depth = m_depth;

    for (amrex::MFIter mfi(levelset); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        auto phi = levelset.array(mfi);

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                phi(i, j, k) = z - (depth - Amp * std::sin(M_PI * x / 2.0));
            });
    }
}

} // namespace amr_wind
