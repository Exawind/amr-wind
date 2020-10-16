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

    pp.query("kappa", m_mo.kappa);
    pp.query("mo_gamma_m", m_mo.gamma_m);
    pp.query("mo_gamma_h", m_mo.gamma_h);
    pp.query("mo_beta_m", m_mo.beta_m);
    pp.query("surface_roughness_z0", m_mo.z0);
    pp.query("normal_direction", m_direction);
    pp.queryarr("gravity", m_gravity);
    AMREX_ASSERT((0 <= m_direction) && (m_direction < AMREX_SPACEDIM));

    if (pp.contains("log_law_height")) {
        m_use_fch = false;
        pp.get("log_law_height", m_mo.zref);
    } else {
        m_use_fch = true;
        amrex::Print()
            << "ABLWallFunction: log_law_height not specified for ABL physics. "
               "Assuming log_law_height = first cell height"
            << std::endl;
    }

    pp.get("reference_temperature", m_mo.ref_temp);

    m_tempflux = true;
    m_mo.surf_temp_flux = 0.0;
    if (pp.contains("surface_temp_flux")) {
        pp.query("surface_temp_flux", m_mo.surf_temp_flux);
    } else if (pp.contains("surface_temp_rate")) {
        m_tempflux = false;
        pp.get("surface_temp_rate", m_surf_temp_rate);
        if (pp.contains("surface_temp_init"))
            pp.get("surface_temp_init", m_surf_temp_init);
        else {
            amrex::Print()
                << "ABLWallFunction: Initial surface temperature not found for "
                   "ABL. Assuming to be equal to the reference temperature "
                << m_mo.ref_temp << std::endl;
            m_surf_temp_init = m_mo.ref_temp;
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

    m_mo.alg_type = m_tempflux ? MOData::HEAT_FLUX : MOData::SURFACE_TEMPERATURE;
    m_mo.gravity = utils::vec_mag(m_gravity.data());
}

void ABLWallFunction::init_log_law_height()
{
    if (m_use_fch) {
        const auto& geom = m_mesh.Geom(0);
        m_mo.zref =
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
                    std::floor((m_mo.zref - first_cell_height) / dz));

    // assuming Z is wall normal direction
    m_ncells_x = dhi.x - dlo.x + 1;
    m_ncells_y = dhi.y - dlo.y + 1;

    amrex::Real zcellN = first_cell_height + (m_z_sample_index)*dz;

    m_coeff_interp[0] = 1.0 - (m_mo.zref - zcellN) / dz;
    m_coeff_interp[1] = 1.0 - m_coeff_interp[0];

    amrex::IntVect lo(AMREX_D_DECL(0, 0, m_z_sample_index));
    amrex::IntVect hi(
        AMREX_D_DECL(m_ncells_x - 1, m_ncells_y - 1, m_z_sample_index));

    m_bx_z_sample.setSmall(lo);
    m_bx_z_sample.setBig(hi);

    // 3 velocity component + potential temperature
    m_store_xy_vel_temp.resize(m_bx_z_sample, 4);
}

void ABLWallFunction::update_umean(
    const VelPlaneAveraging& vpa, const FieldPlaneAveraging& tpa)
{
    const auto& time = m_sim.time();

    if (!m_tempflux)
        m_mo.surf_temp = m_surf_temp_init +
                      m_surf_temp_rate *
                          (time.current_time() - m_surf_temp_rate_tstart) /
                          3600.0;
#if 1
    {
        m_mo.vel_mean[0] = vpa.line_average_interpolated(m_mo.zref, 0);
        m_mo.vel_mean[1] = vpa.line_average_interpolated(m_mo.zref, 1);
        m_mo.vmag_mean = vpa.line_hvelmag_average_interpolated(m_mo.zref);
        m_mo.theta_mean = tpa.line_average_interpolated(m_mo.zref, 0);
    }
    m_mo.update_fluxes();
#else
    computeplanar();
    computeusingheatflux();
#endif
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

        if ((m_mo.zref - zminBox) * (zmaxBox - m_mo.zref) <=
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

    for (int i=0; i < AMREX_SPACEDIM; ++i) m_mo.vel_mean[i] = 0.0;

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

    m_mo.vel_mean[0] = umean0 / numCells;
    m_mo.vel_mean[1] = umean1 / numCells;
    m_mo.vmag_mean = mean_windspd / numCells;
    m_mo.theta_mean = mean_pot_temp / numCells;
}

void ABLWallFunction::computeusingheatflux()
{
    m_mo.update_fluxes();

    auto xy_arr = m_store_xy_vel_temp.array();

    const amrex::Real umean0 = m_mo.vel_mean[0];
    const amrex::Real umean1 = m_mo.vel_mean[1];
    const amrex::Real mean_windspd = m_mo.vmag_mean;
    const amrex::Real mean_pot_temp = m_mo.theta_mean;
    const amrex::Real ref_temp = m_mo.ref_temp;

    const amrex::Real tau_xz = umean0 / m_mo.vmag_mean;
    const amrex::Real tau_yz = umean1 / m_mo.vmag_mean;
    const amrex::Real tau_thetaz = -m_mo.surf_temp_flux;
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
#if 1
    BL_PROFILE("amr-wind::ABLVelWallFunc");
    constexpr bool extrapolate = false;
    constexpr int idim = 2;
    auto& repo = velocity.repo();
    auto& density = repo.get_field("density", rho_state);
    auto& viscosity = repo.get_field("velocity_mueff");
    const int nlevels = repo.num_active_levels();

    amrex::Orientation zlo(amrex::Direction::z, amrex::Orientation::low);

    AMREX_ALWAYS_ASSERT(((velocity.bc_type()[zlo] == BC::wall_model) &&
                         (!repo.mesh().Geom(0).isPeriodic(idim))));

    const amrex::Real c0 = (!extrapolate) ? 1.0 : 1.5;
    const amrex::Real c1 = (!extrapolate) ? 0.0 : -0.5;
    const amrex::Real utau2 = m_wall_func.utau() * m_wall_func.utau();
    const auto& mo = m_wall_func.mo();
    const amrex::Real umeanx = mo.vel_mean[0];
    const amrex::Real umeany = mo.vel_mean[1];
    const amrex::Real wspd_mean = mo.vmag_mean;
    const amrex::Real tau_xz = umeanx / wspd_mean;
    const amrex::Real tau_yz = umeany / wspd_mean;

    for (int lev=0; lev < nlevels; ++lev) {
        const auto& geom = repo.mesh().Geom(lev);
        const auto& domain = geom.Domain();
        amrex::MFItInfo mfi_info{};

        auto& rho_lev = density(lev);
        auto& vold_lev = velocity.state(FieldState::Old)(lev);
        auto& vel_lev = velocity(lev);
        auto& eta_lev = viscosity(lev);

        if (amrex::Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(vel_lev, mfi_info); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.validbox();
            auto varr  = vel_lev.array(mfi);
            auto vold_arr = vold_lev.array(mfi);
            auto den = rho_lev.array(mfi);
            auto eta = eta_lev.array(mfi);

            if (!(bx.smallEnd(idim) == domain.smallEnd(idim))) continue;

            amrex::ParallelFor(
                amrex::bdryLo(bx, idim),
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real mu = c0 * eta(i, j, k) + c1 * eta(i, j, k+1);
                    const amrex::Real uu = vold_arr(i, j, k, 0);
                    const amrex::Real vv = vold_arr(i, j, k, 1);
                    const amrex::Real wspd = std::sqrt(uu * uu + vv * vv);

                    // Dirichlet BC
                    varr(i, j, k-1, 2) = 0.0;

                    // Shear stress BC
                    amrex::Real taux =
                        tau_xz * ((uu - umeanx) * wspd_mean + wspd * umeanx) /
                        (wspd_mean * umeanx);
                    amrex::Real tauy =
                        tau_yz * ((vv - umeany) * wspd_mean + wspd * umeany) /
                        (wspd_mean * umeany);

                    varr(i, j, k-1, 0) = taux * den(i, j, k) * utau2 / mu;
                    varr(i, j, k-1, 1) = tauy * den(i, j, k) * utau2 / mu;
                });
        }
    }

#else
    diffusion::wall_model_bc_moeng(
        velocity, m_wall_func.utau(), rho_state, m_wall_func.instplanar());
#endif
}

ABLTempWallFunc::ABLTempWallFunc(Field&, const ABLWallFunction& wall_fuc)
    : m_wall_func(wall_fuc)
{}

void ABLTempWallFunc::operator()(Field& temperature, const FieldState rho_state)
{
#if 1
    constexpr bool extrapolate = false;
    constexpr int idim = 2;
    auto& repo = temperature.repo();

    // Return early if the user hasn't requested a wall model BC for temperature
    amrex::Orientation zlo(amrex::Direction::z, amrex::Orientation::low);
    if ((temperature.bc_type()[zlo] != BC::wall_model) ||
        repo.mesh().Geom(0).isPeriodic(idim))
        return;

    BL_PROFILE("amr-wind::ABLVelWallFunc");
    auto& velocity = repo.get_field("velocity");
    auto& density = repo.get_field("density", rho_state);
    auto& alpha = repo.get_field("temperature_mueff");
    const int nlevels = repo.num_active_levels();

    const amrex::Real c0 = (!extrapolate) ? 1.0 : 1.5;
    const amrex::Real c1 = (!extrapolate) ? 0.0 : -0.5;
    const auto& mo = m_wall_func.mo();
    const amrex::Real wspd_mean = mo.vmag_mean;
    const amrex::Real theta_mean = mo.theta_mean;
    const amrex::Real theta_surf = mo.surf_temp;
    const amrex::Real term1 = (mo.utau * mo.kappa) / (wspd_mean * mo.phi_h());

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = repo.mesh().Geom(lev);
        const auto& domain = geom.Domain();
        amrex::MFItInfo mfi_info{};

        auto& rho_lev = density(lev);
        auto& vold_lev = velocity.state(FieldState::Old)(lev);
        auto& told_lev = temperature.state(FieldState::Old)(lev);
        auto& theta = temperature(lev);
        auto& eta_lev = alpha(lev);

        if (amrex::Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(theta, mfi_info); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.validbox();
            auto vold_arr = vold_lev.array(mfi);
            auto told_arr = told_lev.array(mfi);
            auto tarr = theta.array(mfi);
            auto den = rho_lev.array(mfi);
            auto eta = eta_lev.array(mfi);

            if (!(bx.smallEnd(idim) == domain.smallEnd(idim))) continue;

            amrex::ParallelFor(
                amrex::bdryLo(bx, idim),
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real alphaT =
                        c0 * eta(i, j, k) + c1 * eta(i, j, k + 1);
                    const amrex::Real uu = vold_arr(i, j, k, 0);
                    const amrex::Real vv = vold_arr(i, j, k, 1);
                    const amrex::Real wspd = std::sqrt(uu * uu + vv * vv);

                    const amrex::Real theta = told_arr(i, j, k);
                    const amrex::Real num1 = (theta - theta_mean) * wspd_mean;
                    const amrex::Real num2 = (theta_mean - theta_surf) * wspd;
                    const amrex::Real tauT = term1 * (num1 + num2);
                    tarr(i, j, k - 1) = den(i, j, k) * tauT / alphaT;
                });
        }
    }
#else
    diffusion::temp_wall_model_bc(
        temperature, m_wall_func.instplanar(), rho_state);
#endif
}

} // namespace amr_wind
