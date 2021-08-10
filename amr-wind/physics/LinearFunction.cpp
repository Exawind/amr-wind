#include "amr-wind/physics/LinearFunction.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/trig_ops.H"

namespace amr_wind {

LinearFunction::LinearFunction(const CFDSim& sim)
    : m_repo(sim.repo()),
      m_velocity(sim.repo().get_field("velocity"))
{
    amrex::ParmParse pp("linearcoefficients");
    pp.query("x", m_xcoef);
    pp.query("y", m_ycoef);
    pp.query("z", m_zcoef);

	sim.repo().declare_nd_field("testscalar", 1, 1, 1);  
}

/** Initialize the velocity and density fields at the beginning of the
 *  simulation.
 */
void LinearFunction::initialize_fields(
    int level, const amrex::Geometry& geom)
{
    using namespace utils;

	auto& testscalar = m_repo.get_field("testscalar");
	auto& velocity = m_velocity(level);

    const auto& problo = geom.ProbLoArray();
    const auto& probhi = geom.ProbHiArray();
    const amrex::Real Lx = probhi[0] - problo[0];
    const amrex::Real Ly = probhi[1] - problo[1];
    const amrex::Real Lz = probhi[2] - problo[2];

    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto& dx = geom.CellSizeArray();
		auto ts = testscalar(level).array(mfi);

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                ts(i, j, k, 0) = x / Lx * m_xcoef + y / Ly * m_ycoef + z / Lz * m_zcoef;
            });
    }
}

} // namespace amr_wind
