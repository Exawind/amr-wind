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
void KOmegaSST<Transport>::parse_model_coeffs()
{

    if ((!this->m_sim.pde_manager().constant_density() ||
         this->m_sim.physics_manager().contains("MultiPhase"))) {
        this->m_include_buoyancy = true;
    }
    this->m_buoyancy_factor = (this->m_include_buoyancy) ? 1.0 : 0.0;

    {
        const std::string coeffs_dict = this->model_name() + "_coeffs";
        amrex::ParmParse pp(coeffs_dict);
        pp.query("beta_star", this->m_beta_star);
        pp.query("alpha1", this->m_alpha1);
        pp.query("alpha2", this->m_alpha2);
        pp.query("beta1", this->m_beta1);
        pp.query("beta2", this->m_beta2);
        pp.query("sigma_k1", this->m_sigma_k1);
        pp.query("sigma_k2", this->m_sigma_k2);
        pp.query("sigma_omega1", this->m_sigma_omega1);
        pp.query("sigma_omega2", this->m_sigma_omega2);
        pp.query("sigma_t", this->m_sigma_t);
    }

    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
    }
}

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
    // cppcheck-suppress constVariable
    auto& repo = mu_turb.repo();

    const int nlevels = repo.num_active_levels();

    auto gradK = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradK, tke);

    auto gradOmega = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradOmega, sdr);

    // This is used for the buoyancy-modified version of the model
    auto gradden = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradden, den);

    const auto& vel = this->m_vel.state(fstate);
    // Compute strain rate into shear production term
    fvm::strainrate(this->m_shear_prod, vel);

    auto& tke_lhs = (this->m_sim).repo().get_field("tke_lhs_src_term");
    tke_lhs.setVal(0.0);
    // cppcheck-suppress constVariable
    auto& sdr_lhs = (this->m_sim).repo().get_field("sdr_lhs_src_term");

    const amrex::Real deltaT = (this->m_sim).time().deltaT();
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        {m_gravity[0], m_gravity[1], m_gravity[2]}};
    const amrex::Real Bfac = this->m_buoyancy_factor;
    const amrex::Real sigmat = this->m_sigma_t;

    for (int lev = 0; lev < nlevels; ++lev) {
        for (amrex::MFIter mfi(mu_turb(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& lam_mu_arr = (*lam_mu)(lev).array(mfi);
            const auto& mu_arr = mu_turb(lev).array(mfi);
            const auto& gradrho_arr = (*gradden)(lev).array(mfi);
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
            const auto& buoy_arr = (this->m_buoy_term(lev)).array(mfi);

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

                    amrex::Real arg2 =
                        amrex::max<amrex::Real>(2.0 * tmp2, tmp3);
                    amrex::Real f2 = std::tanh(arg2 * arg2);

                    mu_arr(i, j, k) =
                        rho_arr(i, j, k) * a1 * tke_arr(i, j, k) /
                        amrex::max(a1 * sdr_arr(i, j, k), tmp4 * f2);

                    // Buoyancy term
                    amrex::Real tmpB =
                        -(gravity[0] * gradrho_arr(i, j, k, 0) +
                          gravity[1] * gradrho_arr(i, j, k, 1) +
                          gravity[2] * gradrho_arr(i, j, k, 2));

                    buoy_arr(i, j, k) = Bfac * tmpB *
                                        (mu_arr(i, j, k) / rho_arr(i, j, k)) /
                                        sigmat;

                    f1_arr(i, j, k) = tmp_f1;

                    diss_arr(i, j, k) = -beta_star * rho_arr(i, j, k) *
                                        tke_arr(i, j, k) * sdr_arr(i, j, k);
                    tke_lhs_arr(i, j, k) = 0.5 * beta_star * rho_arr(i, j, k) *
                                           sdr_arr(i, j, k) * deltaT;

                    shear_prod_arr(i, j, k) = amrex::min<amrex::Real>(
                        mu_arr(i, j, k) * tmp4 * tmp4,
                        10.0 * beta_star * rho_arr(i, j, k) * tke_arr(i, j, k) *
                            sdr_arr(i, j, k));

                    sdr_lhs_arr(i, j, k) = 0.5 * rho_arr(i, j, k) * beta *
                                           sdr_arr(i, j, k) * deltaT;
                    sdr_src_arr(i, j, k) =
                        rho_arr(i, j, k) * alpha * tmp4 * tmp4 +
                        (1.0 - tmp_f1) * 2.0 * rho_arr(i, j, k) * sigma_omega2 *
                            gko / (sdr_arr(i, j, k) + 1e-15);
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
        const auto& repo = deff.repo();
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
        const auto& repo = deff.repo();
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
