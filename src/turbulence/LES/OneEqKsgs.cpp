#include "OneEqKsgs.H"
#include "PDEBase.H"
#include "TurbModelDefs.H"
#include "derive_K.H"
#include "turb_utils.H"
#include "tke/TKE.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace turbulence {

template <typename Transport>
OneEqKsgs<Transport>::OneEqKsgs(CFDSim& sim)
    : TurbModelBase<Transport>(sim), m_vel(sim.repo().get_field("velocity")),
      m_turb_lscale(sim.repo().declare_field("turb_lscale",1, 1, 1)),
      m_shear_prod(sim.repo().declare_field("shear_prod",1, 1, 1)),
      m_buoy_prod(sim.repo().declare_field("buoy_prod",1, 1, 1)),
      m_rho(sim.repo().get_field("density"))
{
    auto& tke_eqn = sim.pde_manager().register_transport_pde(pde::TKE::pde_name());
    m_tke = &(tke_eqn.fields().field);

    //Turbulent length scale field
    this->m_sim.io_manager().register_io_var("turb_lscale");
}

template <typename Transport>
OneEqKsgs<Transport>::~OneEqKsgs() = default;

template <typename Transport>
OneEqKsgsM84<Transport>::OneEqKsgsM84(CFDSim& sim)
    : OneEqKsgs<Transport>(sim),
      m_temperature(sim.repo().get_field("temperature"))
{

    auto& phy_mgr = this->m_sim.physics_manager();
    if (!phy_mgr.contains("ABL")) 
        amrex::Abort("OneEqKsgsM84 model only works with ABL physics");
        
    {
        const std::string coeffs_dict = this->model_name() + "_coeffs";
        amrex::ParmParse pp(coeffs_dict);
        pp.query("Ceps", this->m_Ceps);
        pp.query("Ce", this->m_Ce);
    }

    {
        amrex::ParmParse pp("ABL");
        pp.get("reference_temperature", m_ref_theta);
    }

    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
    }
    
    // TKE source term to be added to PDE
    // turb_utils::inject_turbulence_src_terms(pde::TKE::pde_name(), {"KsgsM84Src"});
}

template <typename Transport>
OneEqKsgsM84<Transport>::~OneEqKsgsM84() = default;

template <typename Transport>
TurbulenceModel::CoeffsDictType OneEqKsgsM84<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{{"Ce", this->m_Ce}, {"Ceps", this->m_Ceps}};
}

template <typename Transport>
void OneEqKsgsM84<Transport>::update_turbulent_viscosity(
    const FieldState fstate)
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::update_turbulent_viscosity")

    auto gradT = this->m_sim.repo().create_scratch_field(AMREX_SPACEDIM,0);
    compute_gradient(*gradT, m_temperature);

    auto& vel = this->m_vel.state(fstate);
    // Compute strain rate into shear production term
    compute_strainrate(this->m_shear_prod, vel);
    
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};
    const amrex::Real beta = 1.0/m_ref_theta;
    
    auto& mu_turb = this->mu_turb();
    const amrex::Real Ce = this->m_Ce;
    auto& den = this->m_rho.state(fstate);
    auto& repo = mu_turb.repo();
    auto& geom_vec = repo.mesh().Geom();

    const int nlevels = repo.num_active_levels();
    for (int lev=0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];

        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real ds = std::cbrt(dx * dy * dz);

        for (amrex::MFIter mfi(mu_turb(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.growntilebox(mu_turb.num_grow());
            const auto& mu_arr = mu_turb(lev).array(mfi);
            const auto& rho_arr = den(lev).const_array(mfi);
            const auto& gradT_arr = (*gradT)(lev).array(mfi);
            const auto& tlscale_arr = (this->m_turb_lscale)(lev).array(mfi);
            const auto& tke_arr = (*this->m_tke)(lev).array(mfi);
            const auto& buoy_prod_arr = (this->m_buoy_prod)(lev).array(mfi);
            const auto& shear_prod_arr = (this->m_shear_prod)(lev).array(mfi);
            
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

                  amrex::Real stratification =
                      (gradT_arr(i,j,k,0) * gravity[0]
                       + gradT_arr(i,j,k,1) * gravity[1]
                       + gradT_arr(i,j,k,2) * gravity[2])*beta;
                  if(stratification > 1e-10)
                      tlscale_arr(i,j,k) =
                          0.76 * stratification * std::sqrt(tke_arr(i,j,k));
                  else
                      tlscale_arr(i,j,k) = ds;
                        
                  mu_arr(i, j, k) =
                      rho_arr(i, j, k) * Ce
                      * tlscale_arr(i,j,k) * std::sqrt(tke_arr(i,j,k));
                  
                  buoy_prod_arr(i,j,k) =
                      2.0 * mu_arr(i,j,k)/(1.0 + 2.0 * tlscale_arr(i,j,k)/ds)
                      * stratification;

                  shear_prod_arr(i,j,k) *= mu_arr(i,j,k);
                  
                    });
        }
    }
    
}

template <typename Transport>
void OneEqKsgsM84<Transport>::update_alphaeff(Field& alphaeff)
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::update_alphaeff")
    
    amrex::Real lam_diff = (this->m_transport).thermal_diffusivity();
    auto& mu_turb = this->m_mu_turb;
    auto& repo = mu_turb.repo();
    auto& geom_vec = repo.mesh().Geom();

    const int nlevels = repo.num_active_levels();
    for (int lev=0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];

        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real ds = std::cbrt(dx * dy * dz);

        for (amrex::MFIter mfi(mu_turb(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.growntilebox(mu_turb.num_grow());
            const auto& muturb_arr = mu_turb(lev).array(mfi);
            const auto& alphaeff_arr = alphaeff(lev).array(mfi);
            const auto& tlscale_arr = this->m_turb_lscale(lev).array(mfi);
            
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        alphaeff_arr(i, j, k) = lam_diff + 
                            muturb_arr(i,j,k)
                            / (1.0 + 2.0 * tlscale_arr(i,j,k)/ds);
                    });
        }
    }
    
}

template <typename Transport>
OneEqKsgsS94<Transport>::OneEqKsgsS94(CFDSim& sim) : OneEqKsgs<Transport>(sim)
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("Ceps", this->m_Ceps);

    // TKE source term to be added to PDE
    // turb_utils::inject_turbulence_src_terms(pde::TKE::pde_name(), {"KsgsS94Src"});
}

template <typename Transport>
OneEqKsgsS94<Transport>::~OneEqKsgsS94() = default;

template <typename Transport>
TurbulenceModel::CoeffsDictType OneEqKsgsS94<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{{"Ceps", this->m_Ceps}};
}

template <typename Transport>
void OneEqKsgsS94<Transport>::update_turbulent_viscosity(
    const FieldState fstate)
{
    BL_PROFILE(
        "amr-wind::" + this->identifier() + "::update_turbulent_viscosity")

    auto& mu_turb = this->mu_turb();
    auto& vel = this->m_vel.state(fstate);
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(OneEqKsgsM84);
INSTANTIATE_TURBULENCE_MODEL(OneEqKsgsS94);

} // namespace amr_wind
