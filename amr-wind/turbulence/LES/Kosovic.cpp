#include <cmath>
#include "amr-wind/turbulence/LES/Kosovic.H"
#include "amr-wind/turbulence/TurbModelDefs.H"
#include "amr-wind/fvm/nonLinearSum.H"
#include "amr-wind/fvm/strainrate.H"
#include "amr-wind/fvm/divergence.H"
#include "AMReX_REAL.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace turbulence {

template <typename Transport>
// cppcheck-suppress uninitMemberVar
Kosovic<Transport>::Kosovic(CFDSim& sim)
    : TurbModelBase<Transport>(sim)
    , m_vel(sim.repo().get_field("velocity"))
    , m_rho(sim.repo().get_field("density"))
    , m_Nij(sim.repo().declare_field("Nij", 9, 1, 1))
    , m_divNij(sim.repo().declare_field("divNij", 3))
{
    amrex::ParmParse pp("Kosovic");
    pp.query("Cb", m_Cb);
    m_Cs = std::sqrt(8 * (1 + m_Cb) / (27 * M_PI * M_PI));
    m_C1 = std::sqrt(960) * m_Cb / (7 * (1 + m_Cb) * m_Sk);
    m_C2 = m_C1;
    pp.query("refMOL", m_refMOL);
    pp.query("surfaceRANS", m_surfaceRANS);
    if (m_surfaceRANS) {
        m_surfaceFactor = 1;
        pp.query("switchLoc", m_switchLoc);
        pp.query("surfaceRANSExp", m_surfaceRANSExp);
    } else {
        m_surfaceFactor = 0;
    }
    pp.query("writeTerms", m_writeTerms);
    if (m_writeTerms) {
        this->m_sim.io_manager().register_io_var("Nij");
        this->m_sim.io_manager().register_io_var("divNij");
    }
    pp.query("LESOff", m_LESTurnOff);
}
template <typename Transport>
void Kosovic<Transport>::update_turbulent_viscosity(
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
    const bool has_terrain =
        this->m_sim.repo().int_field_exists("terrain_blank");
    const auto* m_terrain_blank =
        has_terrain ? &this->m_sim.repo().get_int_field("terrain_blank")
                    : nullptr;
    // Populate strainrate into the turbulent viscosity arrays to avoid creating
    // a temporary buffer
    fvm::strainrate(mu_turb, vel);
    // Non-linear component Nij is computed here and goes into Body Forcing
    fvm::nonlinearsum(m_Nij, vel);
    fvm::divergence(m_divNij, m_Nij);
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];
        const auto& problo = repo.mesh().Geom(lev).ProbLoArray();
        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real ds = std::cbrt(dx * dy * dz);
        const amrex::Real ds_sqr = ds * ds;
        const amrex::Real smag_factor = Cs_sqr * ds_sqr;
        const amrex::Real locMOL = m_refMOL;
        const amrex::Real locLESTurnOff = m_LESTurnOff;
        const amrex::Real locSwitchLoc = m_switchLoc;
        const amrex::Real locSurfaceRANSExp = m_surfaceRANSExp;
        const amrex::Real locSurfaceFactor = m_surfaceFactor;
        const amrex::Real locC1 = m_C1;
        const auto& mu_arrs = mu_turb(lev).arrays();
        const auto& rho_arrs = den(lev).const_arrays();
        const auto& divNij_arrs = (this->m_divNij)(lev).arrays();
        const auto& blank_arrs = has_terrain
                                     ? (*m_terrain_blank)(lev).const_arrays()
                                     : amrex::MultiArray4<const int>();
        amrex::ParallelFor(
            mu_turb(lev),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const amrex::Real rho = rho_arrs[nbx](i, j, k);
                const amrex::Real x3 = problo[2] + (k + 0.5) * dz;
                const amrex::Real fmu = std::exp(-x3 / locSwitchLoc);
                const amrex::Real phiM =
                    (locMOL < 0) ? std::pow(1 - 16 * x3 / locMOL, -0.25)
                                 : 1 + 5 * x3 / locMOL;
                const amrex::Real ransL =
                    std::pow(0.41 * (k + 1) * dz / phiM, 2);
                amrex::Real turnOff = std::exp(-x3 / locLESTurnOff);
                amrex::Real viscosityScale =
                    locSurfaceFactor *
                        (std::pow(1 - fmu, locSurfaceRANSExp) * smag_factor +
                         std::pow(fmu, locSurfaceRANSExp) * ransL) +
                    (1 - locSurfaceFactor) * smag_factor;
                const amrex::Real blankTerrain =
                    (has_terrain) ? 1 - blank_arrs[nbx](i, j, k, 0) : 1.0;
                mu_arrs[nbx](i, j, k) *=
                    rho * viscosityScale * turnOff * blankTerrain;
                amrex::Real stressScale =
                    locSurfaceFactor *
                        (std::pow(1 - fmu, locSurfaceRANSExp) * smag_factor *
                             0.25 * locC1 +
                         std::pow(fmu, locSurfaceRANSExp) * ransL) +
                    (1 - locSurfaceFactor) * smag_factor * 0.25 * locC1;
                divNij_arrs[nbx](i, j, k, 0) *=
                    rho * stressScale * turnOff * blankTerrain;
                divNij_arrs[nbx](i, j, k, 1) *=
                    rho * stressScale * turnOff * blankTerrain;
                divNij_arrs[nbx](i, j, k, 2) *=
                    rho * stressScale * turnOff * blankTerrain;
            });
    }
    amrex::Gpu::streamSynchronize();

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

template <typename Transport>
void Kosovic<Transport>::update_alphaeff(Field& alphaeff)
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::update_alphaeff");

    auto lam_alpha = (this->m_transport).alpha();
    auto& mu_turb = this->m_mu_turb;
    auto& repo = mu_turb.repo();
    const amrex::Real muCoeff = (m_refMOL < 0) ? 3 : 1;
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& muturb_arrs = mu_turb(lev).const_arrays();
        const auto& alphaeff_arrs = alphaeff(lev).arrays();
        const auto& lam_diff_arrs = (*lam_alpha)(lev).const_arrays();
        amrex::ParallelFor(
            mu_turb(lev),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                alphaeff_arrs[nbx](i, j, k) =
                    lam_diff_arrs[nbx](i, j, k) +
                    muCoeff * muturb_arrs[nbx](i, j, k);
            });
    }
    amrex::Gpu::streamSynchronize();

    alphaeff.fillpatch(this->m_sim.time().current_time());
}
template <typename Transport>
void Kosovic<Transport>::parse_model_coeffs()
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("Cs", this->m_Cs);
}

template <typename Transport>
TurbulenceModel::CoeffsDictType Kosovic<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{{"Cb", this->m_Cb}};
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(Kosovic);

} // namespace amr_wind
