#include <AMReX_Orientation.H>

#include "amr-wind/equation_systems/tke/source_terms/KsgsM84Src.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/turbulence/TurbulenceModel.H"

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
    m_CepsGround =  (3.9/0.93) * m_Ceps;

}

KsgsM84Src::~KsgsM84Src() = default;

void KsgsM84Src::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& tlscale_arr = (this->m_turb_lscale)(lev).array(mfi);
    const auto& shear_prod_arr = (this->m_shear_prod)(lev).array(mfi);
    const auto& buoy_prod_arr = (this->m_buoy_prod)(lev).array(mfi);
    const auto& tke_arr = (this->m_tke)(lev).array(mfi);
    const amrex::Real Ceps = this->m_Ceps;
    const amrex::Real CepsGround = this->m_CepsGround;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      src_term(i, j, k) += shear_prod_arr(i,j,k)
          + buoy_prod_arr(i,j,k)
          - Ceps * std::sqrt(tke_arr(i,j,k)) * tke_arr(i,j,k) / tlscale_arr(i,j,k);
    });

    auto& bctype = (this->m_tke).bc_type();
    for (int dir=0; dir < AMREX_SPACEDIM; ++dir) {
      amrex::Orientation olo(dir, amrex::Orientation::low);
      if (bctype[olo] == BC::wall_model) {
        amrex::Box blo = amrex::adjCellLo(bx,dir,0) & bx;
        if (blo.ok()) {
          amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            src_term(i, j, k) += (Ceps - CepsGround) * std::sqrt(tke_arr(i,j,k)) * tke_arr(i,j,k) / tlscale_arr(i,j,k);
          });
        }
      }
    }
    
    
    
}

}
}
}
