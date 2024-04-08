#include <AMReX_GpuContainers.H>
#include <AMReX_GpuQualifiers.H>
#include <cmath>

#include "amr-wind/fvm/gradient.H"
#include "amr-wind/turbulence/LES/AMD.H"
#include "amr-wind/turbulence/TurbModelDefs.H"
#include "amr-wind/utilities/DirectionSelector.H"

#include "AMReX_REAL.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace turbulence {

template <typename Transport>
// cppcheck-suppress uninitMemberVar
AMD<Transport>::AMD(CFDSim& sim)
    : TurbModelBase<Transport>(sim)
    , m_vel(sim.repo().get_field("velocity"))
    , m_temperature(sim.repo().get_field("temperature"))
    , m_rho(sim.repo().get_field("density"))
    , m_pa_temp(m_temperature, sim.time(), m_normal_dir, true)
{
    auto& phy_mgr = this->m_sim.physics_manager();
    if (phy_mgr.contains("ABL")) {
        {
            amrex::ParmParse pp("ABL");
            pp.get("reference_temperature", m_ref_theta);
        }
        {
            amrex::ParmParse pp("incflo");
            pp.queryarr("gravity", m_gravity);
        }
    }
}

template <typename Transport>
void AMD<Transport>::parse_model_coeffs()
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("C_poincare", m_C);
}

template <typename Transport>
void AMD<Transport>::update_turbulent_viscosity(
    const FieldState fstate, const DiffusionType /*unused*/)
{
    BL_PROFILE(
        "amr-wind::" + this->identifier() + "::update_turbulent_viscosity");

    auto& mu_turb = this->mu_turb();
    const auto& repo = mu_turb.repo();
    const auto& vel = m_vel.state(fstate);
    const auto& temp = m_temperature.state(fstate);
    const auto& den = m_rho.state(fstate);
    const auto& geom_vec = repo.mesh().Geom();
    amrex::Real beta = 0.0;
    const int normal_dir = m_normal_dir;
    auto& phy_mgr = this->m_sim.physics_manager();
    if (phy_mgr.contains("ABL")) {
        beta = -m_gravity[normal_dir] / m_ref_theta;
    }

    const amrex::Real C_poincare = m_C;

    auto gradVel = repo.create_scratch_field(AMREX_SPACEDIM * AMREX_SPACEDIM);
    fvm::gradient(*gradVel, vel);
    auto gradT = repo.create_scratch_field(AMREX_SPACEDIM);
    fvm::gradient(*gradT, temp);
    m_pa_temp(); // compute the current plane average
    const auto& tpa_deriv = m_pa_temp.line_deriv();
    amrex::Vector<amrex::Real> tpa_coord(tpa_deriv.size(), 0.0);
    for (int i = 0; i < m_pa_temp.ncell_line(); ++i) {
        tpa_coord[i] = m_pa_temp.xlo() + (0.5 + i) * m_pa_temp.dx();
    }
    amrex::Gpu::DeviceVector<amrex::Real> tpa_deriv_d(tpa_deriv.size(), 0.0);
    amrex::Gpu::DeviceVector<amrex::Real> tpa_coord_d(tpa_coord.size(), 0.0);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, tpa_coord.begin(), tpa_coord.end(),
        tpa_coord_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, tpa_deriv.begin(), tpa_deriv.end(),
        tpa_deriv_d.begin());

    const amrex::Real* p_tpa_coord = tpa_coord_d.data();
    const amrex::Real* p_tpa_deriv = tpa_deriv_d.data();
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];
        const auto& problo = geom.ProbLoArray();
        const amrex::Real nlo = problo[normal_dir];
        const auto& dx = geom.CellSizeArray();

        const auto& gradVel_arrs = (*gradVel)(lev).const_arrays();
        const auto& gradT_arrs = (*gradT)(lev).const_arrays();
        const auto& rho_arrs = den(lev).const_arrays();
        const auto& mu_arrs = mu_turb(lev).arrays();
        const auto& mu_turb_lev = mu_turb(lev);
        amrex::ParallelFor(
            mu_turb_lev,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                auto mu_arr = mu_arrs[nbx];
                const auto rho_arr = rho_arrs[nbx];
                const auto gradVel_arr = gradVel_arrs[nbx];
                const auto gradT_arr = gradT_arrs[nbx];
                mu_arr(i, j, k) =
                    rho_arr(i, j, k) * amd_muvel(
                                           i, j, k, dx, beta, C_poincare,
                                           gradVel_arr, gradT_arr, tpa_deriv_d,
                                           tpa_coord_d, normal_dir, nlo);
            });
    }

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

//! Update the effective thermal diffusivity field
template <typename Transport>
void AMD<Transport>::update_alphaeff(Field& alphaeff)
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::update_alphaeff");

    const auto& repo = alphaeff.repo();
    const auto& geom_vec = repo.mesh().Geom();
    const amrex::Real C_poincare = m_C;
    auto gradVel = repo.create_scratch_field(AMREX_SPACEDIM * AMREX_SPACEDIM);
    fvm::gradient(*gradVel, m_vel);
    auto gradT = repo.create_scratch_field(AMREX_SPACEDIM);
    fvm::gradient(*gradT, m_temperature);

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];

        const auto& dx = geom.CellSizeArray();
        for (amrex::MFIter mfi(alphaeff(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& gradVel_arr = (*gradVel)(lev).const_array(mfi);
            const auto& gradT_arr = (*gradT)(lev).const_array(mfi);
            const auto& rho_arr = m_rho(lev).const_array(mfi);
            const auto& alpha_arr = alphaeff(lev).array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real rho = rho_arr(i, j, k);
                    alpha_arr(i, j, k) = rho * amd_thermal_diff(
                                                   i, j, k, dx, C_poincare,
                                                   gradVel_arr, gradT_arr);
                });
        }
    }
}

template <typename Transport>
TurbulenceModel::CoeffsDictType AMD<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{{"C_poincare", m_C}};
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(AMD);

} // namespace amr_wind
