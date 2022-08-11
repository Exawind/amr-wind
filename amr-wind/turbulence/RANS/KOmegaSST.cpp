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
    pp.query("walldist_type", this->m_walldist_type);
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
    const std::string walldist_type = this->m_walldist_type;

    auto lam_mu = (this->m_transport).mu();
    const auto& den = this->m_rho.state(fstate);
    auto& tke = (*this->m_tke).state(fstate);
    auto& sdr = (*this->m_sdr).state(fstate);
    // cppcheck-suppress constVariable
    auto& repo = mu_turb.repo();

    const auto& geom_vec = repo.mesh().Geom();
    //const auto& dx = geom_vec.CellSizeArray();

    const bool mesh_mapping = this->m_sim.has_mesh_mapping();
    amr_wind::Field const* mesh_fac =
        mesh_mapping ? &((this->m_sim.repo())
                             .get_mesh_mapping_field(amr_wind::FieldLoc::CELL))
                     : nullptr;

    const int nlevels = repo.num_active_levels();

    if (tke.in_uniform_space() && mesh_mapping) {
        tke.to_stretched_space();
    }

    auto gradK = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradK, tke);

    if (sdr.in_uniform_space() && mesh_mapping) {
        sdr.to_stretched_space();
    }

    auto gradOmega = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradOmega, sdr);

    auto& vel = this->m_vel.state(fstate);
    
    // ensure velocity is in stretched mesh space
    if (vel.in_uniform_space() && mesh_mapping) {
        vel.to_stretched_space();
    }

    // Compute strain rate into shear production term
    if (mesh_mapping) {
        shear_prod_to_uniform_space(this->m_shear_prod, vel);
    } else {
        fvm::strainrate(this->m_shear_prod, vel);
    }

    auto& tke_lhs = (this->m_sim).repo().get_field("tke_lhs_src_term");
    tke_lhs.setVal(0.0);
    // cppcheck-suppress constVariable
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
            const auto& vel_arr = vel(lev).array(mfi);
            const auto& wd_arr = (this->m_walldist)(lev).array(mfi);
            const auto& shear_prod_arr = (this->m_shear_prod)(lev).array(mfi);
            const auto& diss_arr = (this->m_diss)(lev).array(mfi);
            const auto& sdr_src_arr = (this->m_sdr_src)(lev).array(mfi);
            const auto& sdr_diss_arr = (this->m_sdr_diss)(lev).array(mfi);
            const auto& f1_arr = (this->m_f1)(lev).array(mfi);
            const auto& tke_lhs_arr = tke_lhs(lev).array(mfi);
            const auto& sdr_lhs_arr = sdr_lhs(lev).array(mfi);
            amrex::Array4<amrex::Real const> fac =
                mesh_mapping ? ((*mesh_fac)(lev).const_array(mfi))
                             : amrex::Array4<amrex::Real const>();

            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const auto& geom = geom_vec[lev];

                    amrex::Real fac_x = mesh_mapping ? fac(i, j, k, 0) : 1.0;
                    amrex::Real fac_y = mesh_mapping ? fac(i, j, k, 1) : 1.0;
                    amrex::Real fac_z = mesh_mapping ? fac(i, j, k, 2) : 1.0;

                    // TODO: Make this general...right now it's z-height only 
                    // assuming that prob-lo is at 0.0
                    amrex::Real zcoord = (k+0.5)*geom.CellSize()[2];
                    amrex::Real zcoord_mapped = this->m_sim.mesh_mapping()->interp_unif_to_nonunif(zcoord, 2);
                    
                    if(walldist_type == "abl"){
                        wd_arr(i, j, k) = mesh_mapping ? zcoord_mapped : zcoord;
                    }

                    if(walldist_type == "channel"){
                        wd_arr(i, j, k) = mesh_mapping ? (zcoord_mapped > 1.0 ? 2.0-zcoord_mapped : zcoord_mapped ) : zcoord;
                    }

                    amrex::Real gko =
                        (gradK_arr(i, j, k, 0) / fac_x *
                             gradOmega_arr(i, j, k, 0) / fac_x +
                         gradK_arr(i, j, k, 1) / fac_y *
                             gradOmega_arr(i, j, k, 1) / fac_y +
                         gradK_arr(i, j, k, 2) / fac_z *
                             gradOmega_arr(i, j, k, 2) / fac_z);

                    amrex::Real cross = 2.0 * rho_arr(i, j, k) * sigma_omega2 * gko /
                                   (sdr_arr(i, j, k) + 1e-15);

                    amrex::Real cdkomega = amrex::max(cross, 1e-10);

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
                        (1.0 - tmp_f1) * cross;

                    sdr_diss_arr(i, j, k) = -rho_arr(i, j, k) * beta *
                                            sdr_arr(i, j, k) * sdr_arr(i, j, k);
                    
                    if(i==16 and j==16 and k<4){
                       amrex::Print()<<shear_prod_arr(i, j, k)<<"\t"<<tmp4<<"\t"<<mu_arr(i,j,k)<<"\t"<<wd_arr(i,j,k)<<"\t"<<vel_arr(i,j,k)<<"\t"<<gko<<std::endl;
                    }
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

template <typename Transport>
void KOmegaSST<Transport>::shear_prod_to_uniform_space(
    Field& shear_prod, const Field& vel)
{
    // Modify shear production while accounting for mesh mapping
    auto grad_vel =
        (this->m_sim).repo().create_scratch_field(9, 0, FieldLoc::CELL);

    fvm::gradient(*grad_vel, vel);

    auto const& mesh_fac =
        (this->m_sim).repo().get_mesh_mapping_field(amr_wind::FieldLoc::CELL);

    for (int lev = 0; lev < (this->m_sim).repo().num_active_levels(); ++lev) {
        for (amrex::MFIter mfi(shear_prod(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& grad_v = (*grad_vel)(lev).array(mfi);
            const auto& fac = (mesh_fac)(lev).array(mfi);
            const auto& str_rt = (shear_prod)(lev).array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    amrex::Real ux = grad_v(i, j, k, 0);
                    amrex::Real uy = grad_v(i, j, k, 1);
                    amrex::Real uz = grad_v(i, j, k, 2);
                    amrex::Real vx = grad_v(i, j, k, 3);
                    amrex::Real vy = grad_v(i, j, k, 4);
                    amrex::Real vz = grad_v(i, j, k, 5);
                    amrex::Real wx = grad_v(i, j, k, 6);
                    amrex::Real wy = grad_v(i, j, k, 7);
                    amrex::Real wz = grad_v(i, j, k, 8);

                    str_rt(i, j, k) = std::sqrt(
                        2.0 * std::pow(ux / fac(i, j, k, 0), 2) +
                        2.0 * std::pow(vy / fac(i, j, k, 1), 2) +
                        2.0 * std::pow(wz / fac(i, j, k, 2), 2) +
                        std::pow(
                            uy / fac(i, j, k, 1) + vx / fac(i, j, k, 0), 2) +
                        std::pow(
                            vz / fac(i, j, k, 2) + wy / fac(i, j, k, 1), 2) +
                        std::pow(
                            wx / fac(i, j, k, 0) + uz / fac(i, j, k, 2), 2));
                });
        }
    }
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(KOmegaSST);

} // namespace amr_wind
