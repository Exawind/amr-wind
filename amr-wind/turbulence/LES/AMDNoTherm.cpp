#include <AMReX_Config.H>
#include <cmath>

#include "amr-wind/fvm/gradient.H"
#include "amr-wind/turbulence/LES/AMDNoTherm.H"
#include "amr-wind/turbulence/TurbModelDefs.H"

#include "AMReX_REAL.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace turbulence {

template <typename Transport>
AMDNoTherm<Transport>::AMDNoTherm(CFDSim& sim)
    : TurbModelBase<Transport>(sim)
    , m_vel(sim.repo().get_field("velocity"))
    , m_rho(sim.repo().get_field("density"))
{}

template <typename Transport>
void AMDNoTherm<Transport>::parse_model_coeffs()
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("C_poincare", m_C);
}

template <typename Transport>
void AMDNoTherm<Transport>::update_turbulent_viscosity(
    const FieldState fstate, const DiffusionType /*unused*/)
{
    BL_PROFILE(
        "amr-wind::" + this->identifier() + "::update_turbulent_viscosity");

    auto& mu_turb = this->mu_turb();
    auto& repo = mu_turb.repo();
    const auto& vel = m_vel.state(fstate);
    const auto& den = m_rho.state(fstate);
    auto gradVel = repo.create_scratch_field(AMREX_SPACEDIM * AMREX_SPACEDIM);
    fvm::gradient(*gradVel, vel);
    const auto& geom_vec = repo.mesh().Geom();

    const amrex::Real C_poincare = this->m_C;

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];

        const auto& dx = geom.CellSizeArray();
        for (amrex::MFIter mfi(mu_turb(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& gradVel_arr = (*gradVel)(lev).array(mfi);
            const auto& mu_arr = mu_turb(lev).array(mfi);
            const auto& rho_arr = den(lev).const_array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real rho = rho_arr(i, j, k);
                    mu_arr(i, j, k) =
                        rho *
                        amd_base_muvel(i, j, k, dx, C_poincare, gradVel_arr);
                });
        }
    }

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

template <typename Transport>
TurbulenceModel::CoeffsDictType AMDNoTherm<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{{"C_poincare", this->m_C}};
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(AMDNoTherm);

} // namespace amr_wind
