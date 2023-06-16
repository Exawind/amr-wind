#include "amr-wind/turbulence/RANS/KOmegaSSTIDDES.H"
#include "amr-wind/turbulence/RANS/KOmegaSSTI.H"
#include "amr-wind/equation_systems/PDEBase.H"
#include "amr-wind/turbulence/TurbModelDefs.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/fvm/strainrate.H"
#include "amr-wind/fvm/vorticity.H"
#include "amr-wind/fvm/vorticity_mag.H"
#include "amr-wind/turbulence/turb_utils.H"
#include "amr-wind/equation_systems/tke/TKE.H"
#include "amr-wind/equation_systems/sdr/SDR.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace turbulence {

template <typename Transport>
KOmegaSSTIDDES<Transport>::~KOmegaSSTIDDES() = default;

template <typename Transport>
KOmegaSSTIDDES<Transport>::KOmegaSSTIDDES(CFDSim& sim)
    : KOmegaSST<Transport>(sim)
    , m_rans_ind(sim.repo().declare_field("rans_indicator", 1, 1, 1))
{}

template <typename Transport>
void KOmegaSSTIDDES<Transport>::parse_model_coeffs()
{
    KOmegaSST<Transport>::parse_model_coeffs();
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("Cdes1", this->m_Cdes1);
    pp.query("Cdes2", this->m_Cdes2);
    pp.query("Cdt1", this->m_Cdt1);
    pp.query("Cdt2", this->m_Cdt2);
    pp.query("Cl", this->m_Cl);
    pp.query("Ct", this->m_Ct);
    pp.query("Cw", this->m_Cw);
    pp.query("kappa", this->m_kappa);
}

template <typename Transport>
TurbulenceModel::CoeffsDictType KOmegaSSTIDDES<Transport>::model_coeffs() const
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
        {"tke_amb", this->m_tke_amb},
        {"sdr_amb", this->m_sdr_amb},
        {"a1", this->m_a1},
        {"Cdes1", this->m_Cdes1},
        {"Cdes2", this->m_Cdes2},
        {"Cdt1", this->m_Cdt1},
        {"Cdt2", this->m_Cdt2},
        {"Cl", this->m_Cl},
        {"Ct", this->m_Ct},
        {"Cw", this->m_Cw},
        {"kappa", this->m_kappa}};
}

template <typename Transport>
void KOmegaSSTIDDES<Transport>::update_turbulent_viscosity(
    const FieldState fstate, const DiffusionType diff_type)
{
    BL_PROFILE(
        "amr-wind::" + this->identifier() + "::update_turbulent_viscosity");

    const amrex::Real beta_star = this->m_beta_star;
    const amrex::Real alpha1 = this->m_alpha1;
    const amrex::Real alpha2 = this->m_alpha2;
    const amrex::Real beta1 = this->m_beta1;
    const amrex::Real beta2 = this->m_beta2;
    const amrex::Real sigma_omega2 = this->m_sigma_omega2;
    const amrex::Real a1 = this->m_a1;
    const amrex::Real tke_amb = this->m_tke_amb;
    const amrex::Real sdr_amb = this->m_sdr_amb;
    const amrex::Real Cdes1 = this->m_Cdes1;
    const amrex::Real Cdes2 = this->m_Cdes2;
    const amrex::Real Cdt1 = this->m_Cdt1;
    const amrex::Real Cdt2 = this->m_Cdt2;
    const amrex::Real Cl = this->m_Cl;
    const amrex::Real Ct = this->m_Ct;
    const amrex::Real Cw = this->m_Cw;
    const amrex::Real kappa = this->m_kappa;

    auto& mu_turb = this->mu_turb();
    auto lam_mu = (this->m_transport).mu();
    const auto& den = this->m_rho.state(fstate);
    const auto& tke = (*this->m_tke).state(fstate);
    const auto& sdr = (*this->m_sdr).state(fstate);
    // cppcheck-suppress constVariable
    auto& repo = mu_turb.repo();
    const auto& geom_vec = repo.mesh().Geom();
    auto& tke_lhs = (this->m_sim).repo().get_field("tke_lhs_src_term");
    tke_lhs.setVal(0.0);
    // cppcheck-suppress constVariable
    auto& sdr_lhs = (this->m_sim).repo().get_field("sdr_lhs_src_term");

    auto gradK = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradK, tke);

    auto gradOmega = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradOmega, sdr);

    const auto& vel = this->m_vel.state(fstate);
    // Compute strain rate into shear production term
    fvm::strainrate(this->m_shear_prod, vel);

    auto vortmag = (this->m_sim.repo()).create_scratch_field(1, 0);
    fvm::vorticity_mag(*vortmag, vel);

    const amrex::Real deltaT = (this->m_sim).time().deltaT();

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];
        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real hmax =
            amrex::max<amrex::Real>(amrex::max<amrex::Real>(dx, dy), dz);

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
            const auto& vortmag_arr = (*vortmag)(lev).const_array(mfi);
            const auto& rans_ind_arr = (this->m_rans_ind)(lev).array(mfi);
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

                    amrex::Real cdkomega = amrex::max<amrex::Real>(
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
                    amrex::Real tmp5 = vortmag_arr(i, j, k);

                    amrex::Real arg1 = amrex::min<amrex::Real>(
                        amrex::max<amrex::Real>(tmp2, tmp3), tmp1);
                    amrex::Real tmp_f1 = std::tanh(arg1 * arg1 * arg1 * arg1);

                    amrex::Real alpha = tmp_f1 * (alpha1 - alpha2) + alpha2;
                    amrex::Real beta = tmp_f1 * (beta1 - beta2) + beta2;

                    amrex::Real arg2 =
                        amrex::max<amrex::Real>(2.0 * tmp2, tmp3);
                    amrex::Real f2 = std::tanh(arg2 * arg2);

                    mu_arr(i, j, k) = rho_arr(i, j, k) * a1 * tke_arr(i, j, k) /
                                      amrex::max<amrex::Real>(
                                          a1 * sdr_arr(i, j, k), tmp4 * f2);

                    f1_arr(i, j, k) = tmp_f1;

                    const amrex::Real alpha_des = 0.25 - wd_arr(i, j, k) / hmax;
                    const amrex::Real denom =
                        (rho_arr(i, j, k) * kappa * kappa * wd_arr(i, j, k) *
                         wd_arr(i, j, k) *
                         std::sqrt(0.5 * (tmp4 * tmp4 + tmp5 * tmp5)));
                    const amrex::Real rdl = lam_mu_arr(i, j, k) / denom;
                    const amrex::Real rdt = mu_arr(i, j, k) / denom;
                    const amrex::Real fl =
                        std::tanh(std::pow(Cl * Cl * rdl, 10));
                    const amrex::Real ft =
                        std::tanh(std::pow(Ct * Ct * rdt, 3));
                    const amrex::Real fe1 =
                        (alpha < 0) ? 2.0 * std::exp(-9.0 * alpha * alpha)
                                    : 2.0 * std::exp(-11.09 * alpha * alpha);
                    const amrex::Real fe2 =
                        1.0 - amrex::max<amrex::Real>(ft, fl);
                    const amrex::Real fe =
                        fe2 * amrex::max<amrex::Real>((fe1 - 1.0), 0.0);
                    const amrex::Real fb = amrex::min<amrex::Real>(
                        2.0 * std::exp(-9.0 * alpha_des * alpha_des), 1.0);
                    const amrex::Real fdt =
                        1.0 - std::tanh(std::pow(Cdt1 * rdt, Cdt2));
                    const amrex::Real fdtilde =
                        amrex::max<amrex::Real>((1.0 - fdt), fb);
                    const amrex::Real cdes = tmp_f1 * (Cdes1 - Cdes2) + Cdes2;
                    const amrex::Real l_les =
                        cdes *
                        amrex::min<amrex::Real>(
                            Cw * amrex::max<amrex::Real>(wd_arr(i, j, k), hmax),
                            hmax);
                    const amrex::Real l_rans = std::sqrt(tke_arr(i, j, k)) /
                                               (beta_star * sdr_arr(i, j, k));
                    const amrex::Real rans_ind = fdtilde * (1.0 + fe);
                    rans_ind_arr(i, j, k) = rans_ind;
                    const amrex::Real l_iddes = amrex::max<amrex::Real>(
                        1.0e-16, rans_ind * l_rans + (1.0 - fdtilde) * l_les);

                    const amrex::Real sqrt_tke_amb = std::sqrt(tke_amb);
                    const amrex::Real l_sst_amb =
                        sqrt_tke_amb /
                        (beta_star * amrex::max<amrex::Real>(1.0e-16, sdr_amb));
                    const amrex::Real l_iddes_amb = amrex::max<amrex::Real>(
                        1.0e-16,
                        rans_ind * l_sst_amb + (1.0 - fdtilde) * l_les);
                    const amrex::Real diss_amb =
                        rho_arr(i, j, k) * tke_amb * sqrt_tke_amb / l_iddes_amb;

                    // For TKE equation:
                    shear_prod_arr(i, j, k) = amrex::min<amrex::Real>(
                        amrex::max<amrex::Real>(
                            mu_arr(i, j, k) * tmp4 * tmp4, 0.0),
                        10.0 * rho_arr(i, j, k) * tke_arr(i, j, k) *
                            std::sqrt(tke_arr(i, j, k)) / l_iddes);

                    diss_arr(i, j, k) = -rho_arr(i, j, k) *
                                            std::sqrt(tke_arr(i, j, k)) *
                                            tke_arr(i, j, k) / l_iddes +
                                        diss_amb;

                    tke_lhs_arr(i, j, k) = 0.5 * rho_arr(i, j, k) *
                                           std::sqrt(tke_arr(i, j, k)) /
                                           l_iddes * deltaT;

                    // For SDR equation:
                    amrex::Real production_omega =
                        rho_arr(i, j, k) * alpha *
                        amrex::min<amrex::Real>(
                            tmp4 * tmp4,
                            10.0 *
                                amrex::max<amrex::Real>(sdr_arr(i, j, k), 0.0) *
                                std::sqrt(tke_arr(i, j, k)) / l_iddes);

                    amrex::Real cross_diffusion =
                        (1.0 - tmp_f1) * 2.0 * rho_arr(i, j, k) * sigma_omega2 *
                        gko / (sdr_arr(i, j, k) + 1e-15);

                    const amrex::Real sdr_diss_amb =
                        beta * rho_arr(i, j, k) * sdr_amb * sdr_amb;

                    if (diff_type == DiffusionType::Crank_Nicolson) {

                        sdr_src_arr(i, j, k) = production_omega;

                        sdr_diss_arr(i, j, k) = cross_diffusion;

                        sdr_lhs_arr(i, j, k) =
                            (rho_arr(i, j, k) * beta * sdr_arr(i, j, k) +
                             0.5 * std::abs(cross_diffusion) /
                                 (sdr_arr(i, j, k) + 1e-15)) *
                            deltaT;

                    } else {
                        sdr_src_arr(i, j, k) =
                            production_omega + cross_diffusion;

                        sdr_diss_arr(i, j, k) = -rho_arr(i, j, k) * beta *
                                                    sdr_arr(i, j, k) *
                                                    sdr_arr(i, j, k) +
                                                sdr_diss_amb;

                        sdr_lhs_arr(i, j, k) = 0.5 * rho_arr(i, j, k) * beta *
                                               sdr_arr(i, j, k) * deltaT;
                    }
                });
        }
    }

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(KOmegaSSTIDDES);

} // namespace amr_wind
