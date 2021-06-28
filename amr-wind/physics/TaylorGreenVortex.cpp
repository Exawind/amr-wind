#include "amr-wind/physics/TaylorGreenVortex.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/trig_ops.H"

namespace amr_wind {

TaylorGreenVortex::TaylorGreenVortex(const CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
    , m_mesh_fac(sim.repo().get_field("mesh_scaling_factor_cc"))
{
    {
        amrex::ParmParse pp("incflo");
        pp.query("density", m_rho);
    }

    // determine probhi based on if mesh is mapped
    {
        amrex::ParmParse pp("geometry");
        if (pp.contains("prob_hi_unmapped")) {
            pp.getarr("prob_hi_unmapped", m_probhi_unmapped);
        }
        else {
            for (int d = 0; d <= AMREX_SPACEDIM; ++d) {
                m_probhi_unmapped[d] = sim.repo().mesh().Geom(0).ProbHiArray()[d];
            }
        }
    }
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
    const amrex::Real Lx = m_probhi_unmapped[0] - problo[0];
    const amrex::Real Ly = m_probhi_unmapped[1] - problo[1];
    const amrex::Real Lz = m_probhi_unmapped[2] - problo[2];

    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto& dx = geom.CellSizeArray();
        auto vel = velocity.array(mfi);
        amrex::Array4<amrex::Real const> const& fac = m_mesh_fac(level).const_array(mfi);

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0]*fac(i, j, k, 0);
                const amrex::Real y = problo[1] + (j + 0.5) * dx[1]*fac(i, j, k, 1);
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2]*fac(i, j, k, 2);

                vel(i, j, k, 0) = std::sin(two_pi() * x / Lx) *
                                  std::cos(two_pi() * y / Ly) *
                                  cos(two_pi() * z / Lz);
                vel(i, j, k, 1) = -std::cos(two_pi() * x / Lx) *
                                  std::sin(two_pi() * y / Ly) *
                                  cos(two_pi() * z / Lz);
                vel(i, j, k, 2) = 0.0;
            });
    }

    m_velocity.is_mesh_mapped() = true;
}

} // namespace amr_wind
