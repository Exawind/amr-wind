#include "amr-wind/turbulence/RANS/KOmegaSST.H"
#include "amr-wind/turbulence/RANS/KOmegaSSTI.H"
#include "amr-wind/equation_systems/PDEBase.H"
#include "amr-wind/turbulence/TurbModelDefs.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/fvm/strainrate.H"
#include "amr-wind/turbulence/turb_utils.H"
#include "amr-wind/equation_systems/tke/TKE.H"
#include "amr-wind/equation_systems/sdr/SDR.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace turbulence {

template <typename Transport>
TurbulenceModel::CoeffsDictType KOmegaSST<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{
        {"beta_star", this->m_beta_star},
        {"alpha1", this->m_alpha1},
        {"alpha2", this->m_alpha2},
        {"beta1", this->m_beta1},
        {"beta2", this->m_beta2},
        {"sigma_k1", this->m_sigma_k1},
        {"sigma_k2", this->m_sigma_k2},
        {"sigma_omega1", this->m_sigma_omega1},
        {"sigma_omega2", this->m_sigma_omega2},
        {"a1", this->m_a1}};
}

template <typename Transport>
void KOmegaSST<Transport>::update_turbulent_viscosity(const FieldState fstate)
{
    BL_PROFILE(
        "amr-wind::" + this->identifier() + "::update_turbulent_viscosity");

    auto& mu_turb = this->mu_turb();
    const amrex::Real beta_star = this->m_beta_star;
    const amrex::Real alpha1 = this->m_alpha1;
    const amrex::Real alpha2 = this->m_alpha2;
    const amrex::Real beta1 = this->m_beta1;
    const amrex::Real beta2 = this->m_beta2;
    const amrex::Real sigma_omega2 = this->m_sigma_omega2;
    const amrex::Real a1 = this->m_a1;

    auto lam_mu = (this->m_transport).mu();
    const auto& den = this->m_rho.state(fstate);
    const auto& tke = (*this->m_tke).state(fstate);
    const auto& sdr = (*this->m_sdr).state(fstate);
    auto& repo = mu_turb.repo();

    const int nlevels = repo.num_active_levels();

    auto gradK = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradK, tke);

    auto gradOmega = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradOmega, sdr);

    const auto& vel = this->m_vel.state(fstate);
    // Compute strain rate into shear production term
    fvm::strainrate(this->m_shear_prod, vel);

    auto& tke_lhs = (this->m_sim).repo().get_field("tke_lhs_src_term");
    tke_lhs.setVal(0.0);
    auto& sdr_lhs = (this->m_sim).repo().get_field("sdr_lhs_src_term");

    const amrex::Real deltaT = (this->m_sim).time().deltaT();

    for (int lev = 0; lev < nlevels; ++lev) {
        for (amrex::MFIter mfi(mu_turb(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& lam_mu_arr = (*lam_mu)(lev).array(mfi);
            const auto& mu_arr = mu_turb(lev).array(mfi);
            const auto& rho_arr = den(lev).const_array(mfi);
            const auto& gradK_arr = (*gradK)(lev).array(mfi);
            const auto& gradOmega_arr = (*gradOmega)(lev).array(mfi);
            const auto& tke_arr = tke(lev).array(mfi);
            const auto& sdr_arr = sdr(lev).array(mfi);
            const auto& wd_arr = (this->m_walldist)(lev).array(mfi);
            const auto& shear_prod_arr = (this->m_shear_prod)(lev).array(mfi);
            const auto& diss_arr = (this->m_diss)(lev).array(mfi);
            const auto& sdr_src_arr = (this->m_sdr_src)(lev).array(mfi);
            const auto& sdr_diss_arr = (this->m_sdr_diss)(lev).array(mfi);
            const auto& f1_arr = (this->m_f1)(lev).array(mfi);
            const auto& tke_lhs_arr = tke_lhs(lev).array(mfi);
            const auto& sdr_lhs_arr = sdr_lhs(lev).array(mfi);

            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    amrex::Real gko =
                        (gradK_arr(i, j, k, 0) * gradOmega_arr(i, j, k, 0) +
                         gradK_arr(i, j, k, 1) * gradOmega_arr(i, j, k, 1) +
                         gradK_arr(i, j, k, 2) * gradOmega_arr(i, j, k, 2));

                    amrex::Real cdkomega = amrex::max(
                        1e-10, 2.0 * rho_arr(i, j, k) * sigma_omega2 * gko /
                                   (sdr_arr(i, j, k) + 1e-15));

                    amrex::Real tmp1 =
                        4.0 * rho_arr(i, j, k) * sigma_omega2 *
                        tke_arr(i, j, k) /
                        (cdkomega * wd_arr(i, j, k) * wd_arr(i, j, k));
                    amrex::Real tmp2 =
                        std::sqrt(tke_arr(i, j, k)) /
                        (beta_star * sdr_arr(i, j, k) * wd_arr(i, j, k) +
                         1e-15);
                    amrex::Real tmp3 =
                        500.0 * lam_mu_arr(i, j, k) /
                        (wd_arr(i, j, k) * wd_arr(i, j, k) * sdr_arr(i, j, k) *
                             rho_arr(i, j, k) +
                         1e-15);
                    amrex::Real tmp4 = shear_prod_arr(i, j, k);

                    amrex::Real arg1 = amrex::min(amrex::max(tmp2, tmp3), tmp1);
                    amrex::Real tmp_f1 = std::tanh(arg1 * arg1 * arg1 * arg1);

                    amrex::Real alpha = tmp_f1 * (alpha1 - alpha2) + alpha2;
                    amrex::Real beta = tmp_f1 * (beta1 - beta2) + beta2;

                    amrex::Real arg2 = amrex::max(2.0 * tmp2, tmp3);
                    amrex::Real f2 = std::tanh(arg2 * arg2);

                    mu_arr(i, j, k) =
                        rho_arr(i, j, k) * a1 * tke_arr(i, j, k) /
                        amrex::max(a1 * sdr_arr(i, j, k), tmp4 * f2);

                    f1_arr(i, j, k) = tmp_f1;

                    diss_arr(i, j, k) = -beta_star * rho_arr(i, j, k) *
                                        tke_arr(i, j, k) * sdr_arr(i, j, k);
                    tke_lhs_arr(i, j, k) = 0.5 * beta_star * rho_arr(i, j, k) *
                                           sdr_arr(i, j, k) * deltaT;

                    shear_prod_arr(i, j, k) = amrex::min(
                        mu_arr(i, j, k) * tmp4 * tmp4,
                        10.0 * beta_star * rho_arr(i, j, k) * tke_arr(i, j, k) *
                            sdr_arr(i, j, k));

                    sdr_lhs_arr(i, j, k) = 0.5 * rho_arr(i, j, k) * beta *
                                           sdr_arr(i, j, k) * deltaT;
                    sdr_src_arr(i, j, k) =
                        rho_arr(i, j, k) * alpha * shear_prod_arr(i, j, k) /
                            amrex::max(mu_arr(i, j, k), 1.0e-16) +
                        (1.0 - tmp_f1) * cdkomega;
                    sdr_diss_arr(i, j, k) = -rho_arr(i, j, k) * beta *
                                            sdr_arr(i, j, k) * sdr_arr(i, j, k);
                });
        }
    }

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

template <typename Transport>
void KOmegaSST<Transport>::update_scalar_diff(
    Field& deff, const std::string& name)
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::update_scalar_diff");

    auto lam_mu = (this->m_transport).mu();
    const auto& mu_turb = this->mu_turb();

    if (name == pde::TKE::var_name()) {
        const amrex::Real sigma_k1 = this->m_sigma_k1;
        const amrex::Real sigma_k2 = this->m_sigma_k2;
        auto& repo = deff.repo();
        const int nlevels = repo.num_active_levels();
        for (int lev = 0; lev < nlevels; ++lev) {
            for (amrex::MFIter mfi(deff(lev)); mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                const auto& lam_mu_arr = (*lam_mu)(lev).array(mfi);
                const auto& mu_arr = mu_turb(lev).array(mfi);
                const auto& f1_arr = (this->m_f1)(lev).array(mfi);
                const auto& deff_arr = deff(lev).array(mfi);
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        deff_arr(i, j, k) =
                            lam_mu_arr(i, j, k) +
                            (f1_arr(i, j, k) * (sigma_k1 - sigma_k2) +
                             sigma_k2) *
                                mu_arr(i, j, k);
                    });
            }
        }

    } else if (name == pde::SDR::var_name()) {
        const amrex::Real sigma_omega1 = this->m_sigma_omega1;
        const amrex::Real sigma_omega2 = this->m_sigma_omega2;
        auto& repo = deff.repo();
        const int nlevels = repo.num_active_levels();
        for (int lev = 0; lev < nlevels; ++lev) {
            for (amrex::MFIter mfi(deff(lev)); mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                const auto& lam_mu_arr = (*lam_mu)(lev).array(mfi);
                const auto& mu_arr = mu_turb(lev).array(mfi);
                const auto& f1_arr = (this->m_f1)(lev).array(mfi);
                const auto& deff_arr = deff(lev).array(mfi);
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        deff_arr(i, j, k) =
                            lam_mu_arr(i, j, k) +
                            (f1_arr(i, j, k) * (sigma_omega1 - sigma_omega2) +
                             sigma_omega2) *
                                mu_arr(i, j, k);
                    });
            }
        }
    } else {
        amrex::Abort(
            "KOmegaSST:update_scalar_diff not implemented for field " + name);
    }
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(KOmegaSST);

} // namespace amr_wind
