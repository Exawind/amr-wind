#include <AMReX_Orientation.H>

#include "amr-wind/equation_systems/tke/source_terms/TKESrc.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/turbulence/TurbulenceModel.H"

namespace amr_wind {
namespace pde {
namespace tke {

TKESrc::TKESrc(const CFDSim& sim)
  : m_shear_prod(sim.repo().get_field("shear_prod")),
    m_diss(sim.repo().get_field("dissipation"))
{}

TKESrc::~TKESrc() = default;

void TKESrc::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& shear_prod_arr = (this->m_shear_prod)(lev).array(mfi);
    const auto& diss_arr = (this->m_diss)(lev).array(mfi);
    
    amrex::Real factor = 1.0;
    if (fstate == FieldState::NPH)
      factor = 0.5;
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
       src_term(i, j, k) += shear_prod_arr(i,j,k) + factor * diss_arr(i,j,k);
    });
    
}

}
}
}
