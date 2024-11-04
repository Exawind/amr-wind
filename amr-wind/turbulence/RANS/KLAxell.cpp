#include "amr-wind/turbulence/RANS/KLAxell.H"
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
KLAxell<Transport>::KLAxell(CFDSim& sim)
    : TurbModelBase<Transport>(sim)
    , m_vel(sim.repo().get_field("velocity"))
    , m_turb_lscale(sim.repo().declare_field("turb_lscale", 1))
    , m_shear_prod(sim.repo().declare_field("shear_prod", 1))
    , m_buoy_prod(sim.repo().declare_field("buoy_prod", 1))
    , m_dissip(sim.repo().declare_field("dissipation", 1))
    , m_rho(sim.repo().get_field("density"))
    , m_temperature(sim.repo().get_field("temperature"))
{
    auto& tke_eqn =
        sim.pde_manager().register_transport_pde(pde::TKE::pde_name());
    m_tke = &(tke_eqn.fields().field);
    auto& phy_mgr = this->m_sim.physics_manager();
    if (!phy_mgr.contains("ABL")) {
        amrex::Abort("KLAxell model only works with ABL physics");
    }
    {
        amrex::ParmParse pp("ABL");
        pp.get("reference_temperature", m_ref_theta);
        pp.get("surface_temp_flux", m_surf_flux);
        pp.query("length_scale_switch", m_lengthscale_switch);
    }

    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
    }

    // TKE source term to be added to PDE
    turb_utils::inject_turbulence_src_terms(
        pde::TKE::pde_name(), {"KransAxell"});
}

template <typename Transport>
void KLAxell<Transport>::parse_model_coeffs()
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("Cmu", this->m_Cmu);
    pp.query("Cmu_prime", this->m_Cmu_prime);
    pp.query("Cb_stable", this->m_Cb_stable);
    pp.query("Cb_unstable", this->m_Cb_unstable);
    pp.query("prandtl", this->m_prandtl);
}

template <typename Transport>
TurbulenceModel::CoeffsDictType KLAxell<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{
        {"Cmu", this->m_Cmu},
        {"Cmu_prime", this->m_Cmu_prime},
        {"Cb_stable", this->m_Cb_stable},
        {"Cb_unstable", this->m_Cb_unstable},
        {"prandtl", this->m_prandtl}};
}

template <typename Transport>
void KLAxell<Transport>::update_turbulent_viscosity(
    const FieldState fstate, const DiffusionType /*unused*/)
{
    BL_PROFILE(
        "amr-wind::" + this->identifier() + "::update_turbulent_viscosity");

    auto gradT = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradT, m_temperature.state(fstate));

    const auto& vel = this->m_vel.state(fstate);
    fvm::strainrate(this->m_shear_prod, vel);

    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};
    const amrex::Real beta = 1.0 / m_ref_theta;
    const amrex::Real Cmu = m_Cmu;
    const amrex::Real Cb_stable = m_Cb_stable;
    const amrex::Real Cb_unstable = m_Cb_unstable;
    auto& mu_turb = this->mu_turb();
    const auto& den = this->m_rho.state(fstate);
    const auto& repo = mu_turb.repo();
    const auto& geom_vec = repo.mesh().Geom();
    const int nlevels = repo.num_active_levels();
    const amrex::Real Rtc = -1.0;
    const amrex::Real Rtmin = -3.0;
    const amrex::Real lambda = 30.0;
    const amrex::Real kappa = 0.41;
    const amrex::Real surf_flux = m_surf_flux;
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real lengthscale_switch = m_lengthscale_switch;
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];
        const auto& problo = repo.mesh().Geom(lev).ProbLoArray();
        const amrex::Real dz = geom.CellSize()[2];
        for (amrex::MFIter mfi(mu_turb(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& mu_arr = mu_turb(lev).array(mfi);
            const auto& rho_arr = den(lev).const_array(mfi);
            const auto& gradT_arr = (*gradT)(lev).array(mfi);
            const auto& tlscale_arr = (this->m_turb_lscale)(lev).array(mfi);
            const auto& tke_arr = (*this->m_tke)(lev).array(mfi);
            const auto& buoy_prod_arr = (this->m_buoy_prod)(lev).array(mfi);
            const auto& shear_prod_arr = (this->m_shear_prod)(lev).array(mfi);
            //! Add terrain components
            const bool has_terrain =
                this->m_sim.repo().int_field_exists("terrain_blank");
            if (has_terrain) {
                const auto* m_terrain_height =
                    &this->m_sim.repo().get_field("terrain_height");
                const auto* m_terrain_blank =
                    &this->m_sim.repo().get_int_field("terrain_blank");
                const auto& ht_arr = (*m_terrain_height)(lev).const_array(mfi);
                const auto& blank_arr =
                    (*m_terrain_blank)(lev).const_array(mfi);
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real stratification =
                            -(gradT_arr(i, j, k, 0) * gravity[0] +
                              gradT_arr(i, j, k, 1) * gravity[1] +
                              gradT_arr(i, j, k, 2) * gravity[2]) *
                            beta;
                        const amrex::Real z = std::max(
                            problo[2] + (k + 0.5) * dz - ht_arr(i, j, k),
                            0.5 * dz);
                        const amrex::Real lscale_s =
                            (lambda * kappa * z) / (lambda + kappa * z);
                        const amrex::Real lscale_b =
                            Cb_stable * std::sqrt(
                                            tke_arr(i, j, k) /
                                            std::max(stratification, tiny));
                        amrex::Real epsilon = std::pow(Cmu, 3) *
                                              std::pow(tke_arr(i, j, k), 1.5) /
                                              (tlscale_arr(i, j, k) + tiny);
                        amrex::Real Rt =
                            std::pow(tke_arr(i, j, k) / epsilon, 2) *
                            stratification;
                        Rt = (Rt > Rtc)
                                 ? Rt
                                 : std::max(
                                       Rt, Rt - std::pow(Rt - Rtc, 2) /
                                                    (Rt + Rtmin - 2 * Rtc));
                        tlscale_arr(i, j, k) =
                            (stratification > 0)
                                ? std::sqrt(
                                      std::pow(lscale_s * lscale_b, 2) /
                                      (std::pow(lscale_s, 2) +
                                       std::pow(lscale_b, 2)))
                                : lscale_s *
                                      std::sqrt(
                                          1.0 -
                                          std::pow(Cmu, 6.0) *
                                              std::pow(Cb_unstable, -2.0) * Rt);
                        tlscale_arr(i, j, k) =
                            (stratification > 0)
                                ? std::min(
                                      tlscale_arr(i, j, k),
                                      std::sqrt(
                                          Cmu * tke_arr(i, j, k) /
                                          stratification))
                                : tlscale_arr(i, j, k);
                        tlscale_arr(i, j, k) = (std::abs(surf_flux) < 1e-5 &&
                                                z <= lengthscale_switch)
                                                   ? lscale_s
                                                   : tlscale_arr(i, j, k);
                        Rt = (std::abs(surf_flux) < 1e-5 &&
                              z <= lengthscale_switch)
                                 ? 0.0
                                 : Rt;
                        const amrex::Real Cmu_Rt =
                            (Cmu + 0.108 * Rt) /
                            (1 + 0.308 * Rt + 0.00837 * std::pow(Rt, 2));
                        mu_arr(i, j, k) = rho_arr(i, j, k) * Cmu_Rt *
                                          tlscale_arr(i, j, k) *
                                          std::sqrt(tke_arr(i, j, k)) *
                                          (1 - blank_arr(i, j, k));
                        const amrex::Real Cmu_prime_Rt = Cmu / (1 + 0.277 * Rt);
                        const amrex::Real muPrime =
                            rho_arr(i, j, k) * Cmu_prime_Rt *
                            tlscale_arr(i, j, k) * std::sqrt(tke_arr(i, j, k)) *
                            (1 - blank_arr(i, j, k));
                        buoy_prod_arr(i, j, k) = -muPrime * stratification;
                        shear_prod_arr(i, j, k) *=
                            shear_prod_arr(i, j, k) * mu_arr(i, j, k);
                    });
            } else {
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real stratification =
                            -(gradT_arr(i, j, k, 0) * gravity[0] +
                              gradT_arr(i, j, k, 1) * gravity[1] +
                              gradT_arr(i, j, k, 2) * gravity[2]) *
                            beta;
                        const amrex::Real z = problo[2] + (k + 0.5) * dz;
                        const amrex::Real lscale_s =
                            (lambda * kappa * z) / (lambda + kappa * z);
                        const amrex::Real lscale_b =
                            Cb_stable * std::sqrt(
                                            tke_arr(i, j, k) /
                                            std::max(stratification, tiny));
                        amrex::Real epsilon = std::pow(Cmu, 3) *
                                              std::pow(tke_arr(i, j, k), 1.5) /
                                              (tlscale_arr(i, j, k) + tiny);
                        amrex::Real Rt =
                            std::pow(tke_arr(i, j, k) / epsilon, 2) *
                            stratification;
                        Rt = (Rt > Rtc)
                                 ? Rt
                                 : std::max(
                                       Rt, Rt - std::pow(Rt - Rtc, 2) /
                                                    (Rt + Rtmin - 2 * Rtc));
                        tlscale_arr(i, j, k) =
                            (stratification > 0)
                                ? std::sqrt(
                                      std::pow(lscale_s * lscale_b, 2) /
                                      (std::pow(lscale_s, 2) +
                                       std::pow(lscale_b, 2)))
                                : lscale_s *
                                      std::sqrt(
                                          1.0 -
                                          std::pow(Cmu, 6.0) *
                                              std::pow(Cb_unstable, -2.0) * Rt);
                        tlscale_arr(i, j, k) =
                            (stratification > 0)
                                ? std::min(
                                      tlscale_arr(i, j, k),
                                      std::sqrt(
                                          Cmu * tke_arr(i, j, k) /
                                          stratification))
                                : tlscale_arr(i, j, k);
                        tlscale_arr(i, j, k) = (std::abs(surf_flux) < 1e-5 &&
                                                z <= lengthscale_switch)
                                                   ? lscale_s
                                                   : tlscale_arr(i, j, k);
                        Rt = (std::abs(surf_flux) < 1e-5 &&
                              z <= lengthscale_switch)
                                 ? 0.0
                                 : Rt;
                        const amrex::Real Cmu_Rt =
                            (Cmu + 0.108 * Rt) /
                            (1 + 0.308 * Rt + 0.00837 * std::pow(Rt, 2));
                        mu_arr(i, j, k) = rho_arr(i, j, k) * Cmu_Rt *
                                          tlscale_arr(i, j, k) *
                                          std::sqrt(tke_arr(i, j, k));
                        const amrex::Real Cmu_prime_Rt = Cmu / (1 + 0.277 * Rt);
                        const amrex::Real muPrime =
                            rho_arr(i, j, k) * Cmu_prime_Rt *
                            tlscale_arr(i, j, k) * std::sqrt(tke_arr(i, j, k));
                        buoy_prod_arr(i, j, k) = -muPrime * stratification;
                        shear_prod_arr(i, j, k) *=
                            shear_prod_arr(i, j, k) * mu_arr(i, j, k);
                    });
            }
        }
    }

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

template <typename Transport>
void KLAxell<Transport>::update_alphaeff(Field& alphaeff)
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::update_alphaeff");
    auto lam_alpha = (this->m_transport).alpha();
    auto& mu_turb = this->m_mu_turb;
    auto& repo = mu_turb.repo();
    auto gradT = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradT, m_temperature);
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};
    const amrex::Real beta = 1.0 / m_ref_theta;
    const amrex::Real Cmu = m_Cmu;
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        for (amrex::MFIter mfi(mu_turb(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& muturb_arr = mu_turb(lev).array(mfi);
            const auto& alphaeff_arr = alphaeff(lev).array(mfi);
            const auto& lam_diff_arr = (*lam_alpha)(lev).array(mfi);
            const auto& tke_arr = (*this->m_tke)(lev).array(mfi);
            const auto& gradT_arr = (*gradT)(lev).array(mfi);
            const auto& tlscale_arr = (this->m_turb_lscale)(lev).array(mfi);
            const amrex::Real Rtc = -1.0;
            const amrex::Real Rtmin = -3.0;
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    amrex::Real stratification =
                        -(gradT_arr(i, j, k, 0) * gravity[0] +
                          gradT_arr(i, j, k, 1) * gravity[1] +
                          gradT_arr(i, j, k, 2) * gravity[2]) *
                        beta;
                    amrex::Real epsilon = std::pow(Cmu, 3) *
                                          std::pow(tke_arr(i, j, k), 1.5) /
                                          tlscale_arr(i, j, k);
                    amrex::Real Rt = std::pow(tke_arr(i, j, k) / epsilon, 2) *
                                     stratification;
                    Rt = (Rt > Rtc) ? Rt
                                    : std::max(
                                          Rt, Rt - std::pow(Rt - Rtc, 2) /
                                                       (Rt + Rtmin - 2 * Rtc));
                    const amrex::Real prandtlRt =
                        (1 + 0.193 * Rt) / (1 + 0.0302 * Rt);
                    alphaeff_arr(i, j, k) =
                        lam_diff_arr(i, j, k) + muturb_arr(i, j, k) / prandtlRt;
                });
        }
    }

    alphaeff.fillpatch(this->m_sim.time().current_time());
}

template <typename Transport>
void KLAxell<Transport>::update_scalar_diff(
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
            "KLAxell:update_scalar_diff not implemented for field " + name);
    }
}

template <typename Transport>
void KLAxell<Transport>::post_advance_work()
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::post_advance_work");
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(KLAxell);

} // namespace amr_wind
