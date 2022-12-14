#include "amr-wind/turbulence/LES/OneEqKsgs.H"
#include "amr-wind/equation_systems/PDEBase.H"
#include "amr-wind/turbulence/TurbModelDefs.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/fvm/strainrate.H"
#include "amr-wind/turbulence/turb_utils.H"
#include "amr-wind/equation_systems/tke/TKE.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace turbulence {

template <typename Transport>
// cppcheck-suppress uninitMemberVar
OneEqKsgs<Transport>::OneEqKsgs(CFDSim& sim)
    : TurbModelBase<Transport>(sim)
    , m_vel(sim.repo().get_field("velocity"))
    , m_turb_lscale(sim.repo().declare_field("turb_lscale", 1, 1, 1))
    , m_shear_prod(sim.repo().declare_field("shear_prod", 1, 1, 1))
    , m_buoy_prod(sim.repo().declare_field("buoy_prod", 1, 1, 1))
    , m_rho(sim.repo().get_field("density"))
{
    auto& tke_eqn =
        sim.pde_manager().register_transport_pde(pde::TKE::pde_name());
    m_tke = &(tke_eqn.fields().field);
}

template <typename Transport>
OneEqKsgs<Transport>::~OneEqKsgs() = default;

template <typename Transport>
OneEqKsgsM84<Transport>::OneEqKsgsM84(CFDSim& sim)
    : OneEqKsgs<Transport>(sim)
    , m_temperature(sim.repo().get_field("temperature"))
{

    auto& phy_mgr = this->m_sim.physics_manager();
    if (!phy_mgr.contains("ABL")) {
        amrex::Abort("OneEqKsgsM84 model only works with ABL physics");
    }

    {
        amrex::ParmParse pp("ABL");
        pp.get("reference_temperature", m_ref_theta);
        pp.query("enable_hybrid_rl_mode", m_hybrid_rl);
    }

    if (m_hybrid_rl) {
        m_sdr = &(sim.repo().declare_field(
            "sdr", 1, (*this->m_tke).num_grow()[0], 1));

        m_sdr->set_default_fillpatch_bc(sim.time());
    }

    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
    }

    // TKE source term to be added to PDE
    turb_utils::inject_turbulence_src_terms(
        pde::TKE::pde_name(), {"KsgsM84Src"});
}

template <typename Transport>
OneEqKsgsM84<Transport>::~OneEqKsgsM84() = default;

template <typename Transport>
void OneEqKsgsM84<Transport>::parse_model_coeffs()
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("Ceps", this->m_Ceps);
    pp.query("Ce", this->m_Ce);
}

template <typename Transport>
TurbulenceModel::CoeffsDictType OneEqKsgsM84<Transport>::model_coeffs() const
{
    // clang-format off
    return TurbulenceModel::CoeffsDictType{
        {"Ce", this->m_Ce}, {"Ceps", this->m_Ceps}};
    // clang-format on
}

template <typename Transport>
void OneEqKsgsM84<Transport>::update_turbulent_viscosity(
    const FieldState fstate)
{
    BL_PROFILE(
        "amr-wind::" + this->identifier() + "::update_turbulent_viscosity");

    auto gradT = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradT, m_temperature.state(fstate));

    const auto& vel = this->m_vel.state(fstate);
    // Compute strain rate into shear production term
    fvm::strainrate(this->m_shear_prod, vel);

    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        {m_gravity[0], m_gravity[1], m_gravity[2]}};
    const amrex::Real beta = 1.0 / m_ref_theta;

    auto& mu_turb = this->mu_turb();
    const amrex::Real Ce = this->m_Ce;
    const auto& den = this->m_rho.state(fstate);
    const auto& repo = mu_turb.repo();
    const auto& geom_vec = repo.mesh().Geom();

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];

        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real ds = std::cbrt(dx * dy * dz);

        for (amrex::MFIter mfi(mu_turb(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
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
                        -(gradT_arr(i, j, k, 0) * gravity[0] +
                          gradT_arr(i, j, k, 1) * gravity[1] +
                          gradT_arr(i, j, k, 2) * gravity[2]) *
                        beta;
                    if (stratification > 1e-10) {
                        tlscale_arr(i, j, k) = amrex::min<amrex::Real>(
                            ds, 0.76 * std::sqrt(
                                           tke_arr(i, j, k) / stratification));
                    } else {
                        tlscale_arr(i, j, k) = ds;
                    }

                    mu_arr(i, j, k) = rho_arr(i, j, k) * Ce *
                                      tlscale_arr(i, j, k) *
                                      std::sqrt(tke_arr(i, j, k));

                    buoy_prod_arr(i, j, k) =
                        -mu_arr(i, j, k) *
                        (1.0 + 2.0 * tlscale_arr(i, j, k) / ds) *
                        stratification;

                    shear_prod_arr(i, j, k) *=
                        shear_prod_arr(i, j, k) * mu_arr(i, j, k);
                });
        }
    }

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

template <typename Transport>
void OneEqKsgsM84<Transport>::update_alphaeff(Field& alphaeff)
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::update_alphaeff");

    auto lam_alpha = (this->m_transport).alpha();
    auto& mu_turb = this->m_mu_turb;
    auto& repo = mu_turb.repo();
    const auto& geom_vec = repo.mesh().Geom();

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];

        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real ds = std::cbrt(dx * dy * dz);

        for (amrex::MFIter mfi(mu_turb(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& muturb_arr = mu_turb(lev).array(mfi);
            const auto& alphaeff_arr = alphaeff(lev).array(mfi);
            const auto& tlscale_arr = this->m_turb_lscale(lev).array(mfi);
            const auto& lam_diff_arr = (*lam_alpha)(lev).array(mfi);

            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    alphaeff_arr(i, j, k) =
                        lam_diff_arr(i, j, k) +
                        muturb_arr(i, j, k) *
                            (1.0 + 2.0 * tlscale_arr(i, j, k) / ds);
                });
        }
    }

    alphaeff.fillpatch(this->m_sim.time().current_time());
}

template <typename Transport>
void OneEqKsgsM84<Transport>::update_scalar_diff(
    Field& deff, const std::string& name)
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::update_scalar_diff");

    if (name == pde::TKE::var_name()) {
        auto& mu_turb = this->mu_turb();
        deff.setVal(0.0);
        field_ops::saxpy(
            deff, 2.0, mu_turb, 0, 0, deff.num_comp(), deff.num_grow());
    } else {
        amrex::Abort(
            "OneEqKsgsM84:update_scalar_diff not implemented for field " +
            name);
    }
}

template <typename Transport>
void OneEqKsgsM84<Transport>::post_advance_work()
{

    if (!m_hybrid_rl) {
        return;
    }

    BL_PROFILE("amr-wind::" + this->identifier() + "::post_advance_work");

    // Update sdr field based on sfs ke

    auto& tke = *(this->m_tke);
    // cppcheck-suppress constVariable
    auto& sdr = *(this->m_sdr);
    const amrex::Real Ce = this->m_Ce;

    auto& repo = tke.repo();
    const auto& geom_vec = repo.mesh().Geom();
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];

        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real ds = std::cbrt(dx * dy * dz);

        for (amrex::MFIter mfi(tke(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.growntilebox();
            const auto& tke_arr = tke(lev).array(mfi);
            const auto& sdr_arr = sdr(lev).array(mfi);

            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    sdr_arr(i, j, k) = std::sqrt(tke_arr(i, j, k)) / (Ce * ds);
                });
        }
    }
}

template <typename Transport>
OneEqKsgsS94<Transport>::OneEqKsgsS94(CFDSim& sim) : OneEqKsgs<Transport>(sim)
{
    // TKE source term to be added to PDE
    turb_utils::inject_turbulence_src_terms(
        pde::TKE::pde_name(), {"KsgsS94Src"});
}

template <typename Transport>
OneEqKsgsS94<Transport>::~OneEqKsgsS94() = default;

template <typename Transport>
void OneEqKsgsS94<Transport>::parse_model_coeffs()
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("Ceps", this->m_Ceps);
}

template <typename Transport>
TurbulenceModel::CoeffsDictType OneEqKsgsS94<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{{"Ceps", this->m_Ceps}};
}

template <typename Transport>
void OneEqKsgsS94<Transport>::update_turbulent_viscosity(
    const FieldState // fstate
    /*unused*/)
{
    BL_PROFILE(
        "amr-wind::" + this->identifier() + "::update_turbulent_viscosity");

    // auto& mu_turb = this->mu_turb();
    // auto& vel = this->m_vel.state(fstate);
}

template <typename Transport>
void OneEqKsgsS94<Transport>::update_scalar_diff(
    Field& deff, const std::string& name)
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::update_scalar_diff");

    if (name == pde::TKE::var_name()) {
        auto& mu_turb = this->mu_turb();
        deff.setVal(0.0);
        field_ops::saxpy(
            deff, 2.0, mu_turb, 0, 0, deff.num_comp(), deff.num_grow());
    } else {
        amrex::Abort(
            "OneEqKsgsM84:update_scalar_diff not implemented for field " +
            name);
    }
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(OneEqKsgsM84);
INSTANTIATE_TURBULENCE_MODEL(OneEqKsgsS94);

} // namespace amr_wind
