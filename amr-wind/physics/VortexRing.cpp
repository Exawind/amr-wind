#include "amr-wind/physics/VortexRing.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/trig_ops.H"
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_FillPatchUtil.H>

namespace amr_wind {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real VorticityTheta::operator()(
    const amrex::Real s,
    const amrex::Real ssq,
    const amrex::Real R,
    const amrex::Real Gamma,
    const amrex::Real alpha) const
{
    if (s <= R) {
		return alpha * Gamma / std::pow(R, 2) * std::exp(- 4 * ssq / (std::pow(R, 2) - ssq));
	} else {
		return 0.0;
	}
}

VortexRing::VortexRing(const CFDSim& sim)
    : m_repo(sim.repo())
	, m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{
	{
    	amrex::ParmParse pp("incflo");
    	pp.query("density", m_rho);
	}

	{
    	amrex::ParmParse pp("vortexring");
    	pp.query("R", m_R);
    	pp.query("Gamma", m_Gamma);
    	pp.query("alpha", m_alpha);
	}

	sim.repo().declare_nd_field("vorticity", 3, 1, 1);    
	sim.repo().declare_nd_field("vectorpotential", 3, 1, 1);    
}

/** Initialize the velocity and density fields at the beginning of the
 *  simulation.
 */
void VortexRing::initialize_fields(
    int level, const amrex::Geometry& geom)
{
    using namespace utils;

    auto& velocity = m_velocity(level);
    auto& density = m_density(level);
	auto& vorticity = m_repo.get_field("vorticity");
    auto& vectorpotential = m_repo.get_field("vectorpotential");

	VorticityTheta vorticity_theta;

    density.setVal(m_rho);

    const auto& problo = geom.ProbLoArray();
    const auto& probhi = geom.ProbHiArray();
    const amrex::Real Lx = probhi[0] - problo[0];
    const amrex::Real Ly = probhi[1] - problo[1];
    const amrex::Real Lz = probhi[2] - problo[2];

	const amrex::Real R = m_R;
	const amrex::Real Gamma = m_Gamma;
	const amrex::Real alpha = m_alpha;

    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& dx = geom.CellSizeArray();
		const auto& nbx = mfi.nodaltilebox();
        auto vort = vorticity(level).array(mfi);

        amrex::ParallelFor(
            nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

                const amrex::Real x = problo[0] + i * dx[0];
                const amrex::Real y = problo[1] + j * dx[1];
                const amrex::Real z = problo[2] + k * dx[2];
                const amrex::Real r = std::sqrt(std::pow(x, 2) + std::pow(y, 2));
                const amrex::Real theta = std::atan2(y, x);
				const amrex::Real ssq = std::pow(z, 2) + std::pow(r - R, 2);
				const amrex::Real s = std::sqrt(ssq);
                vort(i, j, k, 0) = -std::sin(theta) * vorticity_theta(s, ssq, R, Gamma, alpha);
                vort(i, j, k, 1) = std::cos(theta) * vorticity_theta(s, ssq, R, Gamma, alpha);
                vort(i, j, k, 2) = 0.0;
            });
    }

	amrex::LPInfo info;
    auto& mesh = m_velocity.repo().mesh();
    amrex::MLNodeLaplacian linop(mesh.Geom(0,level), mesh.boxArray(0, level), mesh.DistributionMap(0, level), info, {}, 1.0);

    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> bclo;
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> bchi;

    const auto& bctype = m_velocity.bc_type();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (mesh.Geom(0).isPeriodic(dir)) {
            bclo[dir] = amrex::LinOpBCType::Periodic;
            bchi[dir] = amrex::LinOpBCType::Periodic;
        } else {

            switch (bctype[amrex::Orientation(dir, amrex::Orientation::low)]) {
                case BC::pressure_inflow:
                case BC::pressure_outflow: {
                    bclo[dir] = amrex::LinOpBCType::Dirichlet;
                    break;
                }
                default:
                    bclo[dir] = amrex::LinOpBCType::Neumann;
                    break;
            };

            switch (bctype[amrex::Orientation(dir, amrex::Orientation::high)]) {
                case BC::pressure_inflow:
                case BC::pressure_outflow: {
                    bchi[dir] = amrex::LinOpBCType::Dirichlet;
                    break;
                }
                default:
                    bchi[dir] = amrex::LinOpBCType::Neumann;
                    break;
            };
        }
    }

    linop.setDomainBC(bclo,bchi);
    amrex::MLMG mlmg(linop);

    vectorpotential(0).setVal(0.0,0,3,1);

	for(int i=0;i<AMREX_SPACEDIM;++i){
        auto vectorpot = vectorpotential.subview(i,1,level+1);
        auto vort = vorticity.subview(i,1,level+1);
        mlmg.solve(vectorpot.vec_ptrs(), vort.vec_const_ptrs(), 1.0e-6, 0.0);
    }

    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        const auto& dxinv = geom.InvCellSizeArray();
        auto vel = velocity.array(mfi);
        auto psi = vectorpotential(level).array(mfi);
        const amrex::Real facx = amrex::Real(0.25)*dxinv[0];
        const amrex::Real facy = amrex::Real(0.25)*dxinv[1];
        const amrex::Real facz = amrex::Real(0.25)*dxinv[2];

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

                const amrex::Real dpsix_dx = facx * (-psi(i,j,k,0)+psi(i+1,j,k,0)-psi(i,j+1,k,0)+psi(i+1,j+1,k,0)-psi(i,j,k+1,0)+psi(i+1,j,k+1,0)-psi(i,j+1,k+1,0)+psi(i+1,j+1,k+1,0));
                const amrex::Real dpsix_dy = facy * (-psi(i,j,k,0)-psi(i+1,j,k,0)+psi(i,j+1,k,0)+psi(i+1,j+1,k,0)-psi(i,j,k+1,0)-psi(i+1,j,k+1,0)+psi(i,j+1,k+1,0)+psi(i+1,j+1,k+1,0));
                const amrex::Real dpsix_dz = facz * (-psi(i,j,k,0)-psi(i+1,j,k,0)-psi(i,j+1,k,0)-psi(i+1,j+1,k,0)+psi(i,j,k+1,0)+psi(i+1,j,k+1,0)+psi(i,j+1,k+1,0)+psi(i+1,j+1,k+1,0));
                const amrex::Real dpsiy_dx = facx * (-psi(i,j,k,1)+psi(i+1,j,k,1)-psi(i,j+1,k,1)+psi(i+1,j+1,k,1)-psi(i,j,k+1,1)+psi(i+1,j,k+1,1)-psi(i,j+1,k+1,1)+psi(i+1,j+1,k+1,1));
                const amrex::Real dpsiy_dy = facy * (-psi(i,j,k,1)-psi(i+1,j,k,1)+psi(i,j+1,k,1)+psi(i+1,j+1,k,1)-psi(i,j,k+1,1)-psi(i+1,j,k+1,1)+psi(i,j+1,k+1,1)+psi(i+1,j+1,k+1,1));
                const amrex::Real dpsiy_dz = facz * (-psi(i,j,k,1)-psi(i+1,j,k,1)-psi(i,j+1,k,1)-psi(i+1,j+1,k,1)+psi(i,j,k+1,1)+psi(i+1,j,k+1,1)+psi(i,j+1,k+1,1)+psi(i+1,j+1,k+1,1));
                const amrex::Real dpsiz_dx = facx * (-psi(i,j,k,2)+psi(i+1,j,k,2)-psi(i,j+1,k,2)+psi(i+1,j+1,k,2)-psi(i,j,k+1,2)+psi(i+1,j,k+1,2)-psi(i,j+1,k+1,2)+psi(i+1,j+1,k+1,2));
                const amrex::Real dpsiz_dy = facy * (-psi(i,j,k,2)-psi(i+1,j,k,2)+psi(i,j+1,k,2)+psi(i+1,j+1,k,2)-psi(i,j,k+1,2)-psi(i+1,j,k+1,2)+psi(i,j+1,k+1,2)+psi(i+1,j+1,k+1,2));
                const amrex::Real dpsiz_dz = facz * (-psi(i,j,k,2)-psi(i+1,j,k,2)-psi(i,j+1,k,2)-psi(i+1,j+1,k,2)+psi(i,j,k+1,2)+psi(i+1,j,k+1,2)+psi(i,j+1,k+1,2)+psi(i+1,j+1,k+1,2));

                vel(i, j, k, 0) = - (dpsiz_dy - dpsiy_dz);
                vel(i, j, k, 1) = - (dpsix_dz - dpsiz_dx);
                vel(i, j, k, 2) = - (dpsiy_dx - dpsix_dy);

            });
    }
}

} // namespace amr_wind
