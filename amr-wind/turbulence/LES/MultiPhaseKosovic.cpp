#include <cmath>
#include "amr-wind/turbulence/LES/MultiPhaseKosovic.H"
#include "amr-wind/turbulence/TurbModelDefs.H"
#include "amr-wind/fvm/nonLinearSum.H"
#include "amr-wind/fvm/strainrate.H"
#include "amr-wind/fvm/divergence.H"
#include "amr-wind/fvm/gradient.H"
#include "AMReX_REAL.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace turbulence {

template <typename Transport>
// cppcheck-suppress uninitMemberVar
MultiPhaseKosovic<Transport>::MultiPhaseKosovic(CFDSim& sim)
    : TurbModelBase<Transport>(sim)
    , m_vel(sim.repo().get_field("velocity"))
    , m_rho(sim.repo().get_field("density"))
    , m_Nij(sim.repo().declare_field("Nij", 9, 1, 1))
    , m_divNij(sim.repo().declare_field("divNij", 3))
{
    amrex::ParmParse pp("MultiPhaseKosovic");
    pp.query("Cb", m_Cb);
    m_Cs = std::sqrt(8 * (1 + m_Cb) / (27 * M_PI * M_PI));
    m_C1 = std::sqrt(960) * m_Cb / (7 * (1 + m_Cb) * m_Sk);
    m_C2 = m_C1;
    pp.query("writeTerms", m_writeTerms);
    if (m_writeTerms) {
        this->m_sim.io_manager().register_io_var("Nij");
        this->m_sim.io_manager().register_io_var("divNij");
    }
    this->m_sim.io_manager().register_io_var("density_prod");
    {
        amrex::ParmParse pp_incflow("incflo");
        pp_incflow.queryarr("gravity", m_gravity);
    }
    if (!sim.repo().field_exists("vof")) {
        amrex::Abort(
            " Use the single-phase Kosovic model. This model is only for "
            "multiphase");
    }
}
template <typename Transport>
void MultiPhaseKosovic<Transport>::update_turbulent_viscosity(
    const FieldState fstate, const DiffusionType /*unused*/)
{
    BL_PROFILE(
        "amr-wind::" + this->identifier() + "::update_turbulent_viscosity");

    auto& mu_turb = this->mu_turb();
    const auto& repo = mu_turb.repo();
    const auto& vel = m_vel.state(fstate);
    const auto& den = m_rho.state(fstate);
    const auto& geom_vec = repo.mesh().Geom();
    const amrex::Real Cs_sqr = this->m_Cs * this->m_Cs;
    const auto* m_vof = &this->m_sim.repo().get_field("vof");
    // Populate strainrate into the turbulent viscosity arrays to avoid creating
    // a temporary buffer
    fvm::strainrate(mu_turb, vel);
    // Non-linear component Nij is computed here and goes into Body Forcing
    fvm::nonlinearsum(m_Nij, vel);
    fvm::divergence(m_divNij, m_Nij);
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];
        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real ds = std::cbrt(dx * dy * dz);
        const amrex::Real ds_sqr = ds * ds;
        const amrex::Real smag_factor = Cs_sqr * ds_sqr;
        const amrex::Real locC1 = m_C1;
        const auto& mu_arrs = mu_turb(lev).arrays();
        const auto& rho_arrs = den(lev).const_arrays();
        const auto& divNij_arrs = (this->m_divNij)(lev).arrays();
        const auto& vof_arrs = (*m_vof)(lev).const_arrays();
        amrex::ParallelFor(
            mu_turb(lev),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const amrex::Real rho = rho_arrs[nbx](i, j, k);
                //! Turn off turbulence in water as interface not handled with
                //! temperature
                const amrex::Real water_diff =
                    (vof_arrs[nbx](i, j, k) == 0) ? 1 : 0;
                mu_arrs[nbx](i, j, k) *= rho * smag_factor * water_diff;
                //! Non-linear model is not tested in water and interface
                //! effects are approximate Turning off non-linear terms in
                //! water based on density
                amrex::Real water_mixing_length =
                    (vof_arrs[nbx](i, j, k) == 0) ? 1 : 0;
                divNij_arrs[nbx](i, j, k, 0) *=
                    rho * smag_factor * 0.25 * locC1 * water_mixing_length;
                divNij_arrs[nbx](i, j, k, 1) *=
                    rho * smag_factor * 0.25 * locC1 * water_mixing_length;
                divNij_arrs[nbx](i, j, k, 2) *=
                    rho * smag_factor * 0.25 * locC1 * water_mixing_length;
            });
    }
    amrex::Gpu::streamSynchronize();

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

template <typename Transport>
void MultiPhaseKosovic<Transport>::update_alphaeff(Field& alphaeff)
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::update_alphaeff");

    auto lam_alpha = (this->m_transport).alpha();
    auto& mu_turb = this->m_mu_turb;
    auto& repo = mu_turb.repo();
    const amrex::Real muCoeff = 1.0;
    const auto* m_vof = &this->m_sim.repo().get_field("vof");
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& muturb_arrs = mu_turb(lev).const_arrays();
        const auto& alphaeff_arrs = alphaeff(lev).arrays();
        const auto& lam_diff_arrs = (*lam_alpha)(lev).const_arrays();
        const auto& vof_arrs = (*m_vof)(lev).const_arrays();
        amrex::ParallelFor(
            mu_turb(lev),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                //! Current air-water interface setup causes oscillation in
                //! temperature with turbulent diffusion. Turning it off inside
                //! water until a better way to handle it can be found.
                const amrex::Real water_diff =
                    (vof_arrs[nbx](i, j, k) == 0) ? 1 : 0;
                alphaeff_arrs[nbx](i, j, k) =
                    lam_diff_arrs[nbx](i, j, k) +
                    muCoeff * muturb_arrs[nbx](i, j, k) * water_diff;
            });
    }
    amrex::Gpu::streamSynchronize();

    alphaeff.fillpatch(this->m_sim.time().current_time());
}
template <typename Transport>
void MultiPhaseKosovic<Transport>::parse_model_coeffs()
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("Cs", this->m_Cs);
}

template <typename Transport>
TurbulenceModel::CoeffsDictType
MultiPhaseKosovic<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{{"Cb", this->m_Cb}};
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(MultiPhaseKosovic);

} // namespace amr_wind
