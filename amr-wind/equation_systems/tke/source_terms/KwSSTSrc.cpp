#include <AMReX_Orientation.H>

#include "amr-wind/equation_systems/tke/source_terms/KwSSTSrc.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/turbulence/TurbulenceModel.H"

namespace amr_wind {
namespace pde {
namespace tke {

KwSSTSrc::KwSSTSrc(const CFDSim& sim)
    : m_shear_prod(sim.repo().get_field("shear_prod"))
    , m_diss(sim.repo().get_field("dissipation"))
{}

KwSSTSrc::~KwSSTSrc() = default;

void KwSSTSrc::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& shear_prod_arr = (this->m_shear_prod)(lev).array(mfi);
    const auto& diss_arr = (this->m_diss)(lev).array(mfi);

    const amrex::Real factor = (fstate == FieldState::NPH) ? 0.5 : 1.0;
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        src_term(i, j, k) +=
            shear_prod_arr(i, j, k) + factor * diss_arr(i, j, k);
    });
}

} // namespace tke
} // namespace pde
} // namespace amr_wind
