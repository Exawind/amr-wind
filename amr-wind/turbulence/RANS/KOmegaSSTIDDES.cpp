#include "amr-wind/turbulence/RANS/KOmegaSSTIDDES.H"
#include "amr-wind/turbulence/RANS/KOmegaSSTI.H"
#include "amr-wind/equation_systems/PDEBase.H"
#include "amr-wind/turbulence/TurbModelDefs.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/fvm/strainrate.H"
#include "amr-wind/fvm/vorticity.H"
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
{

    {
        const std::string coeffs_dict = this->model_name() + "_coeffs";
        amrex::ParmParse pp(coeffs_dict);
        pp.query("Cdes1", this->m_Cdes1);
        pp.query("Cdes2", this->m_Cdes2);
        pp.query("Cdt1", this->m_Cdt1);
        pp.query("Cdt2", this->m_Cdt2);
        pp.query("Cw", this->m_Cw);
        pp.query("kappa", this->m_kappa);
    }
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
        {"a1", this->m_a1},
        {"Cdes1", this->m_Cdes1},
        {"Cdes2", this->m_Cdes2},
        {"Cdt1", this->m_Cdt1},
        {"Cdt2", this->m_Cdt2},
        {"Cw", this->m_Cw},
        {"kappa", this->m_kappa}};
}

template <typename Transport>
void KOmegaSSTIDDES<Transport>::update_turbulent_viscosity(
    const FieldState fstate)
{
    BL_PROFILE(
        "amr-wind::" + this->identifier() + "::update_turbulent_viscosity");
    (this->m_tke)->fillpatch(this->m_sim.time().current_time());
    (this->m_sdr)->fillpatch(this->m_sim.time().current_time());

    auto gradK = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradK, (this->m_tke)->state(fstate));

    auto gradOmega = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradOmega, (this->m_sdr)->state(fstate));

    auto& vel = this->m_vel.state(fstate);
    // Compute strain rate into shear production term
    fvm::strainrate(this->m_shear_prod, vel);
    // Compute vorticity
    auto vorticity = fvm::vorticity(vel);

    auto& mu_turb = this->mu_turb();
    auto lam_mu = (this->m_transport).mu();
    const amrex::Real beta_star = this->m_beta_star;
    const amrex::Real alpha1 = this->m_alpha1;
    const amrex::Real alpha2 = this->m_alpha2;
    const amrex::Real beta1 = this->m_beta1;
    const amrex::Real beta2 = this->m_beta2;
    const amrex::Real sigma_omega2 = this->m_sigma_omega2;
    const amrex::Real a1 = this->m_a1;
    const amrex::Real Cdes1 = this->m_Cdes1;
    const amrex::Real Cdes2 = this->m_Cdes2;
    const amrex::Real Cdt1 = this->m_Cdt1;
    const amrex::Real Cdt2 = this->m_Cdt2;
    const amrex::Real Cw = this->m_Cw;
    const amrex::Real kappa = this->m_kappa;

    auto& den = this->m_rho.state(fstate);
    auto& repo = mu_turb.repo();
    auto& geom_vec = repo.mesh().Geom();

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];
        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real hmax = amrex::max(amrex::max(dx, dy), dz);

        for (amrex::MFIter mfi(mu_turb(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& lam_mu_arr = (*lam_mu)(lev).array(mfi);
            const auto& mu_arr = mu_turb(lev).array(mfi);
            const auto& rho_arr = den(lev).const_array(mfi);
            const auto& gradK_arr = (*gradK)(lev).array(mfi);
            const auto& gradOmega_arr = (*gradOmega)(lev).array(mfi);
            const auto& tke_arr = (*this->m_tke)(lev).array(mfi);
            const auto& sdr_arr = (*this->m_sdr)(lev).array(mfi);
            const auto& wd_arr = (this->m_walldist)(lev).array(mfi);
            const auto& shear_prod_arr = (this->m_shear_prod)(lev).array(mfi);
            const auto& vort_arr = (*vorticity)(lev).array(mfi);
            const auto& diss_arr = (this->m_diss)(lev).array(mfi);
            const auto& sdr_src_arr = (this->m_sdr_src)(lev).array(mfi);
            const auto& f1_arr = (this->m_f1)(lev).array(mfi);

            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    amrex::Real gko =
                        -(gradK_arr(i, j, k, 0) * gradOmega_arr(i, j, k, 0) +
                          gradK_arr(i, j, k, 1) * gradOmega_arr(i, j, k, 1) +
                          gradK_arr(i, j, k, 2) * gradOmega_arr(i, j, k, 2));

                    amrex::Real cdkomega = amrex::max(
                        1e-10, 2.0 * rho_arr(i, j, k) * sigma_omega2 * gko /
                                   amrex::max(sdr_arr(i, j, k), 1e-15));
                    amrex::Real tmp1 = 4.0 * rho_arr(i, j, k) * sigma_omega2 *
                                       tke_arr(i, j, k) /
                                       (cdkomega * wd_arr(i, j, k));
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
                    amrex::Real tmp5 = vort_arr(i, j, k);

                    amrex::Real tmp_f1 =
                        std::tanh(amrex::min(amrex::max(tmp2, tmp3), tmp1));

                    amrex::Real alpha = tmp_f1 * (alpha1 - alpha2) + alpha2;
                    amrex::Real beta = tmp_f1 * (beta1 - beta2) + beta2;

                    amrex::Real f2 =
                        std::tanh(std::pow(amrex::max(2.0 * tmp2, tmp3), 2));
                    mu_arr(i, j, k) =
                        rho_arr(i, j, k) * a1 * tke_arr(i, j, k) /
                        amrex::max(a1 * sdr_arr(i, j, k), tmp4 * f2);

                    f1_arr(i, j, k) = tmp_f1;

                    const amrex::Real alpha_des = 0.25 - wd_arr(i, j, k) / hmax;
                    const amrex::Real fb = amrex::min(
                        2.0 * std::exp(-9.0 * alpha_des * alpha_des), 1.0);
                    const amrex::Real rdt =
                        (mu_arr(i, j, k) /
                         (rho_arr(i, j, k) * kappa * kappa * wd_arr(i, j, k) *
                          wd_arr(i, j, k) *
                          std::sqrt(0.5 * (tmp4 * tmp4 + tmp5 * tmp5))));
                    const amrex::Real fdt =
                        1.0 - std::tanh(std::pow(Cdt1 * rdt, Cdt2));
                    const amrex::Real fdtilde = amrex::max((1.0 - fdt), fb);
                    const amrex::Real cdes = tmp_f1 * (Cdes1 - Cdes2) + Cdes2;
                    const amrex::Real l_les =
                        cdes *
                        amrex::min(Cw * amrex::max(wd_arr(i, j, k), hmax), hmax);
                    const amrex::Real l_rans = std::sqrt(tke_arr(i, j, k)) /
                                               (beta_star * sdr_arr(i, j, k));
                    const amrex::Real l_iddes =
                        fdtilde * (l_rans - l_les) + l_les;
                    diss_arr(i, j, k) = -std::sqrt(tke_arr(i, j, k)) *
                                        tke_arr(i, j, k) / l_iddes;

                    sdr_src_arr(i, j, k) =
                        rho_arr(i, j, k) *
                        (alpha * tmp4 * tmp4 -
                         beta * sdr_arr(i, j, k) * sdr_arr(i, j, k) +
                         2.0 * (1 - tmp_f1) * sigma_omega2 * gko /
                             (sdr_arr(i, j, k) + 1e-15));

                    shear_prod_arr(i, j, k) = amrex::max(
                        mu_arr(i, j, k) * tmp4 * tmp4,
                        10.0 * beta_star * rho_arr(i, j, k) * tke_arr(i, j, k) *
                            sdr_arr(i, j, k));
                });
        }
    }

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(KOmegaSSTIDDES);

} // namespace amr_wind
