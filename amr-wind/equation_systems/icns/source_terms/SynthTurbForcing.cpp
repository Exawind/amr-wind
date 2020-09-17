#include "amr-wind/equation_systems/icns/source_terms/SynthTurbForcing.H"
#include "amr-wind/utilities/PlaneAveraging.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/trig_ops.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"

namespace amr_wind {
namespace pde {
namespace icns {

/** Synthetic Turbulence forcing term
 *
 *
 *
 */
SynthTurbForcing::SynthTurbForcing(const CFDSim& sim)
    : m_turb_force(sim.repo().get_field("synth_turb_forcing"))
{
//TODO: Check that SyntheticTurbulence physics is enabled
}


SynthTurbForcing::~SynthTurbForcing() = default;

void SynthTurbForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& turb_force_arr = m_turb_force(lev).array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        src_term(i, j, k, 0) += turb_force_arr(i,j,k,0);
        src_term(i, j, k, 1) += turb_force_arr(i,j,k,1);
        src_term(i, j, k, 2) += turb_force_arr(i,j,k,2);
    });
}

} // namespace icns
} // namespace pde
} // namespace amr_wind
