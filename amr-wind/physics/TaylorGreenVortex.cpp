#include "amr-wind/physics/TaylorGreenVortex.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/trig_ops.H"

namespace amr_wind {

TaylorGreenVortex::TaylorGreenVortex(const CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{
    amrex::ParmParse pp("incflo");
    pp.query("density", m_rho);
}

/** Initialize the velocity and density fields at the beginning of the
 *  simulation.
 */
void TaylorGreenVortex::initialize_fields(
    int level, const amrex::Geometry& geom)
{
    using namespace utils;

    auto& velocity = m_velocity(level);
    auto& density = m_density(level);

    density.setVal(m_rho);

    const auto& problo = geom.ProbLoArray();
    const auto& probhi = geom.ProbHiArray();
    const amrex::Real Lx = probhi[0] - problo[0];
    const amrex::Real Ly = probhi[1] - problo[1];
    const amrex::Real Lz = probhi[2] - problo[2];

    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto& dx = geom.CellSizeArray();
        auto vel = velocity.array(mfi);

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                vel(i, j, k, 0) = std::sin(two_pi() * x / Lx) *
                                  std::cos(two_pi() * y / Ly) *
                                  cos(two_pi() * z / Lz);
                vel(i, j, k, 1) = -std::cos(two_pi() * x / Lx) *
                                  std::sin(two_pi() * y / Ly) *
                                  cos(two_pi() * z / Lz);
                vel(i, j, k, 2) = 0.0;
            });
    }
}

} // namespace amr_wind
