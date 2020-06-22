#include "KsgsM84Src.H"
#include "CFDSim.H"
#include "TurbulenceModel.H"

namespace amr_wind {
namespace pde {
namespace tke {

KsgsM84Src::KsgsM84Src(const CFDSim& sim)
  : m_turb_lscale(sim.repo().get_field("turb_lscale")),
    m_shear_prod(sim.repo().get_field("shear_prod")),
    m_buoy_prod(sim.repo().get_field("buoy_prod")),
    m_tke(sim.repo().get_field("tke"))
{
    AMREX_ALWAYS_ASSERT(sim.turbulence_model().model_name() == "OneEqKsgsM84");
    auto coeffs = sim.turbulence_model().model_coeffs();
    m_Ceps = coeffs["Ceps"];
}

KsgsM84Src::~KsgsM84Src() = default;

void KsgsM84Src::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& tlscale_arr = (this->m_turb_lscale)(lev).array(mfi);
    const auto& shear_prod_arr = (this->m_shear_prod)(lev).array(mfi);
    const auto& buoy_prod_arr = (this->m_buoy_prod)(lev).array(mfi);
    const auto& tke_arr = (this->m_tke)(lev).array(mfi);
    const amrex::Real Ceps = this->m_Ceps;
    
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      src_term(i, j, k) += shear_prod_arr(i,j,k)
          + buoy_prod_arr(i,j,k)
          + Ceps * std::sqrt(tke_arr(i,j,k))
          * tke_arr(i,j,k) / tlscale_arr(i,j,k);
    });
    
    
}

}
}
}
