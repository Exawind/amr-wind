#include "amr-wind/wind_energy/ABLWallFunction.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/diffusion/diffusion.H"

#include <cmath>

#include "AMReX_ParmParse.H"
#include "AMReX_Print.H"
#include "AMReX_ParallelDescriptor.H"

namespace amr_wind {

ABLWallFunction::ABLWallFunction(const CFDSim& sim)
    : m_sim(sim), m_mesh(sim.mesh())
{
    amrex::ParmParse pp("ABL");

    pp.query("kappa", m_kappa);
    pp.query("mo_gamma_m", m_gamma_m);
    pp.query("mo_gamma_h", m_gamma_h);
    pp.query("mo_beta_m", m_beta_m);
    pp.query("surface_roughness_z0", m_z0);
    pp.query("normal_direction", m_direction);
    pp.queryarr("gravity", m_gravity);
    AMREX_ASSERT((0 <= m_direction) && (m_direction < AMREX_SPACEDIM));

    if (pp.contains("log_law_height")) {
        m_use_fch = false;
        pp.get("log_law_height", m_log_law_height);
    } else {
        m_use_fch = true;
        amrex::Print()
            << "ABLWallFunction: log_law_height not specified for ABL physics. "
               "Assuming log_law_height = first cell height"
            << std::endl;
    }

    pp.get("reference_temperature", m_ref_temp);

    m_tempflux = true;
    m_surf_temp_flux = 0.0;
    if (pp.contains("surface_temp_flux")) {
        pp.query("surface_temp_flux", m_surf_temp_flux);
    } else if (pp.contains("surface_temp_rate")) {
        m_tempflux = false;
        pp.get("surface_temp_rate", m_surf_temp_rate);
        if (pp.contains("surface_temp_init"))
            pp.get("surface_temp_init", m_surf_temp_init);
        else {
            amrex::Print()
                << "ABLWallFunction: Initial surface temperature not found for "
                   "ABL. Assuming to be equal to the reference temperature "
                << m_ref_temp << std::endl;
            m_surf_temp_init = m_ref_temp;
        }
        if (pp.contains("surface_temp_rate_tstart"))
            pp.get("surface_temp_rate_tstart", m_surf_temp_rate_tstart);
        else {
            amrex::Print()
                << "ABLWallFunction: Surface temperature heating/cooling start "
                   "time (surface_temp_rate_tstart) not found for ABL. "
                   "Assuming zero."
                << m_surf_temp_rate_tstart << std::endl;
        }

    } else {
        amrex::Print() << "ABLWallFunction: Neither surface_temp_flux nor "
                          "surface_temp_rate specified for ABL physics. "
                          "Assuming Neutral Stratification"
                       << std::endl;
    }
}

//! Return Monin-Obukhov Similarity function psi_m
amrex::Real ABLWallFunction::mo_psi_m(amrex::Real zeta)
{

    if (zeta > 0) {
        return -m_gamma_m * zeta;
    } else {
        amrex::Real x = std::sqrt(std::sqrt(1 - m_beta_m * zeta));
        return 2.0 * std::log(0.5 * (1.0 + x)) + log(0.5 * (1 + x * x)) -
               2.0 * std::atan(x) + utils::half_pi();
    }
}

//! Return Monin-Obukhov Similarity function psi_h
amrex::Real ABLWallFunction::mo_psi_h(amrex::Real zeta)
{

    if (zeta > 0) {
        return -m_gamma_h * zeta;
    } else {
        amrex::Real x = std::sqrt(1 - m_beta_m * zeta);
        return std::log(0.5 * (1 + x));
    }
}

void ABLWallFunction::init_log_law_height()
{
    if (m_use_fch) {
        const auto& geom = m_mesh.Geom(0);
        m_log_law_height =
            (geom.ProbLo(m_direction) + 0.5 * geom.CellSize(m_direction));
    }

    const auto& geom = m_mesh.Geom();

    amrex::Box const& domain = geom[m_mesh.finestLevel()].Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    const amrex::Real dz = geom[m_mesh.finestLevel()].CellSize(2);
    amrex::Real first_cell_height =
        geom[m_mesh.finestLevel()].ProbLo(2) + 0.5 * dz;
    m_z_sample_index =
        dlo.z + static_cast<int>(
                    std::floor((m_log_law_height - first_cell_height) / dz));

    // assuming Z is wall normal direction
    m_ncells_x = dhi.x - dlo.x + 1;
    m_ncells_y = dhi.y - dlo.y + 1;

    amrex::Real zcellN = first_cell_height + (m_z_sample_index)*dz;

    m_coeff_interp[0] = 1.0 - (m_log_law_height - zcellN) / dz;
    m_coeff_interp[1] = 1.0 - m_coeff_interp[0];

    amrex::IntVect lo(AMREX_D_DECL(0, 0, m_z_sample_index));
    amrex::IntVect hi(
        AMREX_D_DECL(m_ncells_x - 1, m_ncells_y - 1, m_z_sample_index));

    m_bx_z_sample.setSmall(lo);
    m_bx_z_sample.setBig(hi);

    // 3 velocity component + potential temperature
    m_store_xy_vel_temp.resize(m_bx_z_sample, 4);
}

void ABLWallFunction::update_umean()
{
    const auto& time = m_sim.time();

    if (!m_tempflux)
        m_surf_temp = m_surf_temp_init +
                      m_surf_temp_rate *
                          (time.current_time() - m_surf_temp_rate_tstart) /
                          3600.0;

    computeplanar();
    computeusingheatflux();
}

void ABLWallFunction::computeplanar()
{

    const auto& frepo = m_sim.repo();
    static constexpr int idir = 2;
    const int maxlev = m_mesh.finestLevel();
    const auto& geom = m_mesh.Geom(maxlev);
    auto const& problo = geom.ProbLoArray();
    const amrex::Real dz = geom.CellSize(idir);

    auto& velf = frepo.get_field("velocity", amr_wind::FieldState::New);
    auto& tempf = frepo.get_field("temperature", amr_wind::FieldState::New);

    m_store_xy_vel_temp.setVal<amrex::RunOn::Device>(0.0);

    auto xy_arr = m_store_xy_vel_temp.array();

    const amrex::Real coeff_interp0 = m_coeff_interp[0];
    const amrex::Real coeff_interp1 = m_coeff_interp[1];

    for (amrex::MFIter mfi(velf(maxlev)); mfi.isValid(); ++mfi) {
        const auto& bx = mfi.validbox();

        const auto dlo = amrex::lbound(bx);
        const auto dhi = amrex::ubound(bx);

        amrex::Real zminBox = problo[2] + dz * (dlo.z);
        amrex::Real zmaxBox = problo[2] + dz * (dhi.z);

        if ((m_log_law_height - zminBox) * (zmaxBox - m_log_law_height) <=
            0.0) {
            continue;
        }

        auto vel = velf(maxlev).array(mfi);
        auto temp = tempf(maxlev).array(mfi);

        const amrex::IntVect lo(dlo.x, dlo.y, m_z_sample_index);
        const amrex::IntVect hi(dhi.x, dhi.y, m_z_sample_index);
        const amrex::Box z_sample_bx(lo, hi);

        amrex::ParallelFor(
            z_sample_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                xy_arr(i, j, k, 0) = coeff_interp0 * vel(i, j, k, 0) +
                                       coeff_interp1 * vel(i, j, k + 1, 0);
                xy_arr(i, j, k, 1) = coeff_interp0 * vel(i, j, k, 1) +
                                       coeff_interp1 * vel(i, j, k + 1, 1);
                xy_arr(i, j, k, 2) = coeff_interp0 * vel(i, j, k, 2) +
                                       coeff_interp1 * vel(i, j, k + 1, 2);
                xy_arr(i, j, k, 3) = coeff_interp0 * temp(i, j, k) +
                                       coeff_interp1 * temp(i, j, k + 1);
            });
    }

    amrex::Real numCells = static_cast<amrex::Real>(m_ncells_x * m_ncells_y);

    amrex::ParallelDescriptor::ReduceRealSum(
        m_store_xy_vel_temp.dataPtr(), m_ncells_x * m_ncells_y * 4);

    std::fill(m_umean.begin(), m_umean.end(), 0.0);

    amrex::Real umean0 = 0.0;
    amrex::Real umean1 = 0.0;
    amrex::Real mean_windspd = 0.0;
    amrex::Real mean_pot_temp = 0.0;

    amrex::Loop(
        m_bx_z_sample,
        [=, &umean0, &umean1, &mean_windspd, &mean_pot_temp](int i, int j, int k) noexcept {
            umean0 += xy_arr(i, j, k, 0);
            umean1 += xy_arr(i, j, k, 1);
            mean_windspd += std::sqrt(
                xy_arr(i, j, k, 0) * xy_arr(i, j, k, 0) +
                xy_arr(i, j, k, 1) * xy_arr(i, j, k, 1));
            mean_pot_temp += xy_arr(i, j, k, 3);

        });

    m_umean[0] = umean0 / numCells;
    m_umean[1] = umean1 / numCells;
    m_mean_windspeed = mean_windspd / numCells;
    m_mean_pot_temp = mean_pot_temp / numCells;
}

void ABLWallFunction::computeusingheatflux()
{

    amrex::Real g = utils::vec_mag(m_gravity.data());
    amrex::Real zeta = 0.0;
    amrex::Real m_utau_iter = 0.0;

    // Initialize variables
    m_psi_m = 0.0;
    m_psi_h = 0.0;
    m_utau = m_kappa * m_mean_windspeed /
             (std::log(m_log_law_height / m_z0) - m_psi_m);

    int iter = 0;
    do {
        m_utau_iter = m_utau;
        if (m_tempflux) {
            m_surf_temp = m_surf_temp_flux *
                              (std::log(m_log_law_height / m_z0) - m_psi_h) /
                              (m_utau * m_kappa) +
                          m_mean_pot_temp;
        } else {
            m_surf_temp_flux = -(m_mean_pot_temp - m_surf_temp) * m_utau *
                               m_kappa /
                               (std::log(m_log_law_height / m_z0) - m_psi_h);
        }
        m_obukhov_length = -m_utau * m_utau * m_utau * m_mean_pot_temp /
                           (m_kappa * g * m_surf_temp_flux);
        zeta = m_log_law_height / m_obukhov_length;
        m_psi_m = mo_psi_m(zeta);
        m_psi_h = mo_psi_h(zeta);
        m_utau = m_kappa * m_mean_windspeed /
                 (std::log(m_log_law_height / m_z0) - m_psi_m);
        iter += 1;
    } while ((std::abs(m_utau_iter - m_utau) > 1e-5) && iter <= m_max_iter);

    auto xy_arr = m_store_xy_vel_temp.array();

    const amrex::Real umean0 = m_umean[0];
    const amrex::Real umean1 = m_umean[1];
    const amrex::Real mean_windspd = m_mean_windspeed;
    const amrex::Real mean_pot_temp = m_mean_pot_temp;
    const amrex::Real ref_temp = m_ref_temp;

    const amrex::Real tau_xz = umean0 / m_mean_windspeed;
    const amrex::Real tau_yz = umean1 / m_mean_windspeed;
    const amrex::Real tau_thetaz = -m_surf_temp_flux;
    const amrex::Real denom1 = mean_windspd * (mean_pot_temp - ref_temp);

    amrex::ParallelFor(
        m_bx_z_sample, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const amrex::Real inst_wind_speed = std::sqrt(
                xy_arr(i, j, k, 0) * xy_arr(i, j, k, 0) +
                xy_arr(i, j, k, 1) * xy_arr(i, j, k, 1));

            xy_arr(i, j, k, 0) =
                tau_xz *
                ((xy_arr(i, j, k, 0) - umean0) * mean_windspd +
                 inst_wind_speed * umean0) /
                (mean_windspd * umean0);

            xy_arr(i, j, k, 1) =
                tau_yz *
                ((xy_arr(i, j, k, 1) - umean1) * mean_windspd +
                 inst_wind_speed * umean1) /
                (mean_windspd * umean1);

            const amrex::Real num1 =
                (xy_arr(i, j, k, 3) - mean_pot_temp) * mean_windspd;
            const amrex::Real num2 = inst_wind_speed * (mean_pot_temp - ref_temp);

            xy_arr(i, j, k, 3) = tau_thetaz * (num1 + num2) / denom1;
        });
}

ABLVelWallFunc::ABLVelWallFunc(Field&, const ABLWallFunction& wall_func)
    : m_wall_func(wall_func)
{}

void ABLVelWallFunc::operator()(Field& velocity, const FieldState rho_state)
{

    diffusion::wall_model_bc_moeng(
        velocity, m_wall_func.utau(), rho_state, m_wall_func.instplanar());
}

ABLTempWallFunc::ABLTempWallFunc(Field&, const ABLWallFunction& wall_fuc)
    : m_wall_func(wall_fuc)
{}

void ABLTempWallFunc::operator()(Field& temperature, const FieldState rho_state)
{

    diffusion::temp_wall_model_bc(
        temperature, m_wall_func.instplanar(), rho_state);
}

} // namespace amr_wind
