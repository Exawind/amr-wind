#include "TaylorGreenVortex.H"
#include "CFDSim.H"
#include "AMReX_ParmParse.H"
#include "trig_ops.H"

namespace amr_wind {


TaylorGreenVortex::TaylorGreenVortex(const CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{
    amrex::ParmParse pp("incflo");
    pp.query("ro_0", m_rho);
}

/** Initialize the velocity and density fields at the beginning of the
 *  simulation.
 */
void TaylorGreenVortex::initialize_fields(
    int level,
    const amrex::Geometry& geom)
{
    using namespace utils;

    auto& velocity = m_velocity(level);
    auto& density = m_density(level);

    density.setVal(m_rho);
    
    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        const auto& dx = geom.CellSizeArray();
        const auto& problo = geom.ProbLoArray();
        auto vel = velocity.array(mfi);

        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
            const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
            const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
            vel(i,j,k,0) =  std::sin(two_pi()*x) * std::cos(two_pi()*y) * cos(two_pi()*z);
            vel(i,j,k,1) = -std::cos(two_pi()*x) * std::sin(two_pi()*y) * cos(two_pi()*z);
            vel(i,j,k,2) = 0.0;
        });
    }
}

} // namespace amr_wind
