#include "OneEqKsgs.H"
#include "PDEBase.H"
#include "TurbModelDefs.H"
#include "turb_utils.H"
#include "tke/TKE.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace turbulence {

template <typename Transport>
OneEqKsgs<Transport>::OneEqKsgs(CFDSim& sim)
    : TurbModelBase<Transport>(sim), m_vel(sim.repo().get_field("velocity"))
{
    sim.pde_manager().register_transport_pde(pde::TKE::pde_name());
}

template <typename Transport>
OneEqKsgs<Transport>::~OneEqKsgs() = default;

template <typename Transport>
OneEqKsgsM84<Transport>::OneEqKsgsM84(CFDSim& sim) : OneEqKsgs<Transport>(sim)
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("Ceps", this->m_Ceps);

    // TKE source term to be added to PDE
    // turb_utils::inject_turbulence_src_terms(pde::TKE::pde_name(), {"KsgsM84Src"});
}

template <typename Transport>
OneEqKsgsM84<Transport>::~OneEqKsgsM84() = default;

template <typename Transport>
TurbulenceModel::CoeffsDictType OneEqKsgsM84<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{{"Ceps", this->m_Ceps}};
}

template <typename Transport>
void OneEqKsgsM84<Transport>::update_turbulent_viscosity(
    const FieldState fstate)
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::update_turbulent_viscosity")

    auto& mu_turb = this->mu_turb();
    auto& vel = this->m_vel.state(fstate);
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
