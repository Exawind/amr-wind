#include <cmath>

#include "amr-wind/turbulence/LES/Smagorinsky.H"
#include "amr-wind/turbulence/TurbModelDefs.H"
#include "amr-wind/fvm/strainrate.H"
#include "AMReX_REAL.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"

namespace amr_wind {
namespace turbulence {

template <typename Transport>
// cppcheck-suppress uninitMemberVar
Smagorinsky<Transport>::Smagorinsky(CFDSim& sim)
    : TurbModelBase<Transport>(sim)
    , m_vel(sim.repo().get_field("velocity"))
    , m_rho(sim.repo().get_field("density"))
{}

template <typename Transport>
void Smagorinsky<Transport>::update_turbulent_viscosity(
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

    const bool is_vof = this->m_sim.repo().field_exists("vof");
    amrex::Real rho_min = 1e9;

    const auto& vof = this->m_sim.repo().get_field("vof");
    if (is_vof) {
        rho_min = amr_wind::field_ops::global_min_magnitude(m_rho);
    }

    // const auto& vof_arr = (*m_vof)(lev).const_array(mfi);

    // Populate strainrate into the turbulent viscosity arrays to avoid creating
    // a temporary buffer
    fvm::strainrate(mu_turb, vel);

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];

        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real ds = std::cbrt(dx * dy * dz);
        const amrex::Real ds_sqr = ds * ds;
        const amrex::Real smag_factor = Cs_sqr * ds_sqr;

        for (amrex::MFIter mfi(mu_turb(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& mu_arr = mu_turb(lev).array(mfi);
            const auto& rho_arr = den(lev).const_array(mfi);

            const auto& vof_arr = vof(lev).const_array(mfi);

            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    amrex::Real rho;

                    const bool is_air = vof_arr(i, j, k) < 1.;
                    //
                    if (is_vof) {
                        rho = amrex::min(rho_arr(i, j, k), rho_min);
                    } else {
                        rho = rho_arr(i, j, k);
                    }

                    mu_arr(i, j, k) *= rho * smag_factor;
                });
        }
    }

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

template <typename Transport>
void Smagorinsky<Transport>::parse_model_coeffs()
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("Cs", this->m_Cs);
}

template <typename Transport>
TurbulenceModel::CoeffsDictType Smagorinsky<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{{"Cs", this->m_Cs}};
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(Smagorinsky);

} // namespace amr_wind
