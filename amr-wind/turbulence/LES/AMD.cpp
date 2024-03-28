#include <cmath>

#include "amr-wind/fvm/gradient.H"
#include "amr-wind/turbulence/LES/AMD.H"
#include "amr-wind/turbulence/TurbModelDefs.H"

#include "AMReX_REAL.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace turbulence {

template <typename Transport>
// cppcheck-suppress uninitMemberVar
AMD<Transport>::AMD(CFDSim& sim)
    : TurbModelBase<Transport>(sim)
    , m_vel(sim.repo().get_field("velocity"))
    , m_temperature(sim.repo().get_field("temperature"))
    , m_rho(sim.repo().get_field("density"))
{
    auto& phy_mgr = this->m_sim.physics_manager();
    if (phy_mgr.contains("ABL")) {
        {
            amrex::ParmParse pp("ABL");
            pp.get("reference_temperature", m_ref_theta);
        }
        {
            amrex::ParmParse pp("incflo");
            pp.queryarr("gravity", m_gravity);
        }
    }
}

template <typename Transport>
void AMD<Transport>::parse_model_coeffs()
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("C_poincare", m_C);
}

template <typename Transport>
void AMD<Transport>::update_turbulent_viscosity(
    const FieldState fstate, const DiffusionType /*unused*/)
{
    BL_PROFILE(
        "amr-wind::" + this->identifier() + "::update_turbulent_viscosity");

    auto& mu_turb = this->mu_turb();
    const auto& repo = mu_turb.repo();
    const auto& vel = m_vel.state(fstate);
    const auto& temp = m_temperature.state(fstate);
    const auto& den = m_rho.state(fstate);
    const auto& geom_vec = repo.mesh().Geom();
    amrex::Real beta = 0.0;
    auto& phy_mgr = this->m_sim.physics_manager();
    if (phy_mgr.contains("ABL")) {
        beta = -m_gravity[2] / m_ref_theta;
    }

    const amrex::Real C_poincare = this->m_C;

    auto gradVel = repo.create_scratch_field(AMREX_SPACEDIM * AMREX_SPACEDIM);
    fvm::gradient(*gradVel, vel);
    auto gradT = repo.create_scratch_field(AMREX_SPACEDIM);
    fvm::gradient(*gradT, temp);
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];

        const auto& dx = geom.CellSizeArray();
        for (amrex::MFIter mfi(mu_turb(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& gradVel_arr = (*gradVel)(lev).array(mfi);
            const auto& gradT_arr = (*gradT)(lev).array(mfi);
            const auto& mu_arr = mu_turb(lev).array(mfi);
            const auto& rho_arr = den(lev).const_array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real rho = rho_arr(i, j, k);
                    mu_arr(i, j, k) = rho * amd_muvel(
                                                i, j, k, dx, beta, C_poincare,
                                                gradVel_arr, gradT_arr);
                });
        }
    }

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

//! Update the effective thermal diffusivity field
template <typename Transport>
void AMD<Transport>::update_alphaeff(Field& alphaeff)
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::update_alphaeff");

    const auto& repo = alphaeff.repo();
    const auto& geom_vec = repo.mesh().Geom();
    const amrex::Real C_poincare = this->m_C;
    auto gradVel = repo.create_scratch_field(AMREX_SPACEDIM * AMREX_SPACEDIM);
    fvm::gradient(*gradVel, m_vel);
    auto gradT = repo.create_scratch_field(AMREX_SPACEDIM);
    fvm::gradient(*gradT, m_temperature);

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];

        const auto& dx = geom.CellSizeArray();
        for (amrex::MFIter mfi(alphaeff(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& gradVel_arr = (*gradVel)(lev).array(mfi);
            const auto& gradT_arr = (*gradT)(lev).array(mfi);
            const auto& alpha_arr = alphaeff(lev).array(mfi);
            const auto& rho_arr = m_rho(lev).const_array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real rho = rho_arr(i, j, k);
                    alpha_arr(i, j, k) = rho * amd_thermal_diff(i, j, k, dx, C_poincare,gradVel_arr, gradT_arr);
                });
        }
    }

}

template <typename Transport>
TurbulenceModel::CoeffsDictType AMD<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{{"C_poincare", this->m_C}};
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(AMD);

} // namespace amr_wind
