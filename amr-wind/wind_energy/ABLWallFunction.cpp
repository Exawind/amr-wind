#include "amr-wind/wind_energy/ABLWallFunction.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/diffusion/diffusion.H"
#include "amr-wind/wind_energy/ShearStress.H"

#include <cmath>

#include "AMReX_ParmParse.H"
#include "AMReX_Print.H"
#include "AMReX_ParallelDescriptor.H"

namespace amr_wind {

ABLWallFunction::ABLWallFunction(const CFDSim& sim)
    : m_sim(sim), m_mesh(sim.mesh())
{
    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
    }

    amrex::ParmParse pp("ABL");

    pp.query("kappa", m_mo.kappa);
    pp.query("mo_gamma_m", m_mo.gamma_m);
    pp.query("mo_gamma_h", m_mo.gamma_h);
    pp.query("mo_beta_m", m_mo.beta_m);
    pp.query("mo_beta_h", m_mo.beta_h);
    pp.query("surface_roughness_z0", m_mo.z0);
    pp.query("normal_direction", m_direction);
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

    if (pp.contains("surface_temp_flux")) {
        pp.query("surface_temp_flux", m_mo.surf_temp_flux);
    } else if (pp.contains("surface_temp_rate")) {
        m_tempflux = false;
        pp.get("surface_temp_rate", m_surf_temp_rate);
        if (pp.contains("surface_temp_init")) {
            pp.get("surface_temp_init", m_surf_temp_init);
        } else {
            amrex::Print()
                << "ABLWallFunction: Initial surface temperature not found for "
                   "ABL. Assuming to be equal to the reference temperature "
                << m_mo.ref_temp << std::endl;
            m_surf_temp_init = m_mo.ref_temp;
        }
        if (pp.contains("surface_temp_rate_tstart")) {
            pp.get("surface_temp_rate_tstart", m_surf_temp_rate_tstart);
        } else {
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

    if (pp.contains("inflow_outflow_mode")) {
        pp.query("inflow_outflow_mode", m_inflow_outflow);
        if (m_inflow_outflow) {
            pp.query("wf_velocity", m_wf_vel);
            pp.query("wf_vmag", m_wf_vmag);
            pp.query("wf_theta", m_wf_theta);
            amrex::Print() << "ABLWallFunction: Inflow/Outflow mode is turned "
                              "on. Please make sure wall shear stress type is "
                              "set to local."
                           << std::endl;
        }
    }

    m_mo.alg_type = m_tempflux ? MOData::ThetaCalcType::HEAT_FLUX
                               : MOData::ThetaCalcType::SURFACE_TEMPERATURE;
    m_mo.gravity = utils::vec_mag(m_gravity.data());
}

void ABLWallFunction::init_log_law_height()
{
    if (m_use_fch) {
        const auto& geom = m_mesh.Geom(0);
        m_mo.zref = 0.5 * geom.CellSize(m_direction);
    }
}

void ABLWallFunction::update_umean(
    const VelPlaneAveragingFine& vpa, const FieldPlaneAveragingFine& tpa)
{
    const auto& time = m_sim.time();

    if (!m_tempflux) {
        m_mo.surf_temp =
            m_surf_temp_init +
            m_surf_temp_rate *
                amrex::max<amrex::Real>(
                    time.current_time() - m_surf_temp_rate_tstart, 0.0) /
                3600.0;
    }

    if (m_inflow_outflow) {
        m_mo.vel_mean[0] = m_wf_vel[0];
        m_mo.vel_mean[1] = m_wf_vel[1];
        m_mo.vmag_mean = m_wf_vmag;
        m_mo.Su_mean = 0.0; // TODO: need to fill this correctly
        m_mo.Sv_mean = 0.0; // TODO: need to fill this correctly
        m_mo.theta_mean = m_wf_theta;
    } else {
        m_mo.vel_mean[0] = vpa.line_average_interpolated(m_mo.zref, 0);
        m_mo.vel_mean[1] = vpa.line_average_interpolated(m_mo.zref, 1);
        m_mo.vmag_mean = vpa.line_hvelmag_average_interpolated(m_mo.zref);
        m_mo.Su_mean = vpa.line_su_average_interpolated(m_mo.zref);
        m_mo.Sv_mean = vpa.line_sv_average_interpolated(m_mo.zref);
        m_mo.theta_mean = tpa.line_average_interpolated(m_mo.zref, 0);
    }

    m_mo.update_fluxes();
}

void ABLWallFunction::update_tflux(const amrex::Real tflux)
{
    m_mo.surf_temp_flux = tflux;
}

ABLVelWallFunc::ABLVelWallFunc(
    Field& /*unused*/, const ABLWallFunction& wall_func)
    : m_wall_func(wall_func)
{
    amrex::ParmParse pp("ABL");
    pp.query("wall_shear_stress_type", m_wall_shear_stress_type);
    m_wall_shear_stress_type = amrex::toLower(m_wall_shear_stress_type);

    if (m_wall_shear_stress_type == "constant" ||
        m_wall_shear_stress_type == "local" ||
        m_wall_shear_stress_type == "schumann" ||
        m_wall_shear_stress_type == "donelan" ||
        m_wall_shear_stress_type == "moeng") {
        amrex::Print() << "Shear Stress model: " << m_wall_shear_stress_type
                       << std::endl;
    } else {
        amrex::Abort("Shear Stress wall model input mistake");
    }
}

template <typename ShearStress>
void ABLVelWallFunc::wall_model(
    Field& velocity, const FieldState rho_state, const ShearStress& tau)
{
    BL_PROFILE("amr-wind::ABLVelWallFunc");

    constexpr int idim = 2;
    const auto& repo = velocity.repo();
    const auto& temperature = repo.get_field("temperature");
    const auto& density = repo.get_field("density", rho_state);
    const auto& viscosity = repo.get_field("velocity_mueff");
    const int nlevels = repo.num_active_levels();

    amrex::Orientation zlo(amrex::Direction::z, amrex::Orientation::low);
    amrex::Orientation zhi(amrex::Direction::z, amrex::Orientation::high);
    if (velocity.bc_type()[zhi] == BC::wall_model) {
        amrex::Abort("ABL wall models are not applicable to a zhi BC");
    }
    if (velocity.bc_type()[zlo] != BC::wall_model) {
        return;
    }
    // Using multi-roughness only when available in registry
    const bool has_variable_surf_temp = repo.field_exists("surf_temp");
    const auto* m_surf_temp =
        has_variable_surf_temp ? &repo.get_field("surf_temp") : nullptr;
    const bool has_variable_surf_roughness = repo.field_exists("lowerz0");
    const auto* m_terrainz0 =
        has_variable_surf_roughness ? &repo.get_field("lowerz0") : nullptr;
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = repo.mesh().Geom(lev);
        const auto& domain = geom.Domain();
        amrex::MFItInfo mfi_info{};
        const amrex::Real dz = geom.CellSize()[idim];
        const auto& rho_lev = density(lev);
        const auto& vold_lev = velocity.state(FieldState::Old)(lev);
        const auto& temperatureold_lev =
            temperature.state(FieldState::Old)(lev);
        auto& vel_lev = velocity(lev);
        const auto& eta_lev = viscosity(lev);
        if (amrex::Gpu::notInLaunchRegion()) {
            mfi_info.SetDynamic(true);
        }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(vel_lev, mfi_info); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.validbox();
            const auto& varr = vel_lev.array(mfi);
            const auto& vold_arr = vold_lev.const_array(mfi);
            const auto& temperatureold_arr =
                temperatureold_lev.const_array(mfi);
            const auto& den = rho_lev.const_array(mfi);
            const auto& eta = eta_lev.const_array(mfi);
            const auto& z0_arr = (has_variable_surf_roughness)
                                     ? (*m_terrainz0)(lev).const_array(mfi)
                                     : amrex::Array4<double>();
            const auto& surf_temp_arr =
                (has_variable_surf_temp) ? (*m_surf_temp)(lev).const_array(mfi)
                                         : amrex::Array4<double>();
            const amrex::Real louis_bm = m_louis_bm;
            const amrex::Real louis_dm = m_louis_dm;
            if (bx.smallEnd(idim) == domain.smallEnd(idim) &&
                velocity.bc_type()[zlo] == BC::wall_model) {
                amrex::ParallelFor(
                    amrex::bdryLo(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real mu = eta(i, j, k);
                        const amrex::Real uu = vold_arr(i, j, k, 0);
                        const amrex::Real vv = vold_arr(i, j, k, 1);
                        const amrex::Real wspd = std::sqrt(uu * uu + vv * vv);
                        amrex::Real ustar = 1;
                        // Using an explicit stratification correction to avoid
                        // MOL iterations Can be changed with multi-roughness in
                        // future Ref: J. F. Louis "A PARAMETRIC MODEL OF
                        // VERTICAL EDDY FLUXES IN THE ATMOSPHERE"
                        if (has_variable_surf_roughness &&
                            has_variable_surf_temp) {
                            const amrex::Real theta_surface =
                                surf_temp_arr(i, j, k, 1);
                            const amrex::Real theta =
                                temperatureold_arr(i, j, k);
                            const amrex::Real rib =
                                (9.81 * dz * (theta - theta_surface)) /
                                (theta * wspd * wspd);
                            const amrex::Real louis_cm =
                                7.4 * 0.41 * 0.41 /
                                std::pow(
                                    std::log(0.5 * dz / z0_arr(i, j, k)), 2) *
                                louis_bm *
                                std::sqrt(0.5 * dz / z0_arr(i, j, k));
                            const amrex::Real Fm =
                                (rib > 0)
                                    ? 1 / std::pow(1 + louis_dm * rib, 2)
                                    : 1 - louis_bm * rib /
                                              (1 +
                                               louis_cm *
                                                   std::sqrt(std::abs(rib)));
                            ustar = wspd * 0.41 /
                                    std::log(0.5 * dz / z0_arr(i, j, k)) *
                                    std::sqrt(Fm);
                        }
                        // Dirichlet BC
                        varr(i, j, k - 1, 2) = 0.0;

                        // Shear stress BC
                        varr(i, j, k - 1, 0) =
                            (has_variable_surf_roughness)
                                ? std::pow(ustar, 2) * uu / wspd *
                                      den(i, j, k) / mu
                                : tau.calc_vel_x(uu, wspd) * den(i, j, k) / mu;
                        varr(i, j, k - 1, 1) =
                            (has_variable_surf_roughness)
                                ? std::pow(ustar, 2) * vv / wspd *
                                      den(i, j, k) / mu
                                : tau.calc_vel_y(vv, wspd) * den(i, j, k) / mu;
                    });
            }
        }
    }
}

void ABLVelWallFunc::operator()(Field& velocity, const FieldState rho_state)
{
    const auto& mo = m_wall_func.mo();

    if (m_wall_shear_stress_type == "moeng") {

        auto tau = ShearStressMoeng(mo);
        wall_model(velocity, rho_state, tau);

    } else if (m_wall_shear_stress_type == "constant") {

        auto tau = ShearStressConstant(mo);
        wall_model(velocity, rho_state, tau);

    } else if (m_wall_shear_stress_type == "local") {

        auto tau = ShearStressLocal(mo);
        wall_model(velocity, rho_state, tau);

    } else if (m_wall_shear_stress_type == "schumann") {

        auto tau = ShearStressSchumann(mo);
        wall_model(velocity, rho_state, tau);

    } else if (m_wall_shear_stress_type == "donelan") {

        auto tau = ShearStressDonelan(mo);
        wall_model(velocity, rho_state, tau);
    }
}

ABLTempWallFunc::ABLTempWallFunc(
    Field& /*unused*/, const ABLWallFunction& wall_fuc)
    : m_wall_func(wall_fuc)
{
    amrex::ParmParse pp("ABL");
    pp.query("wall_shear_stress_type", m_wall_shear_stress_type);
    m_wall_shear_stress_type = amrex::toLower(m_wall_shear_stress_type);
}

template <typename HeatFlux>
void ABLTempWallFunc::wall_model(
    Field& temperature, const FieldState rho_state, const HeatFlux& tau)
{
    constexpr int idim = 2;
    auto& repo = temperature.repo();

    // Return early if the user hasn't requested a wall model BC for temperature
    amrex::Orientation zlo(amrex::Direction::z, amrex::Orientation::low);
    amrex::Orientation zhi(amrex::Direction::z, amrex::Orientation::high);
    if (temperature.bc_type()[zhi] == BC::wall_model) {
        amrex::Abort("ABL wall models are not applicable to a zhi BC");
    }
    if (temperature.bc_type()[zlo] != BC::wall_model) {
        return;
    }

    BL_PROFILE("amr-wind::ABLTempWallFunc");
    auto& velocity = repo.get_field("velocity");
    const auto& density = repo.get_field("density", rho_state);
    const auto& alpha = repo.get_field("temperature_mueff");
    const int nlevels = repo.num_active_levels();
    const bool has_variable_surf_temp = repo.field_exists("surf_temp");
    const auto* m_surf_temp =
        has_variable_surf_temp ? &repo.get_field("surf_temp") : nullptr;
    const bool has_variable_surf_roughness = repo.field_exists("lowerz0");
    const auto* m_terrainz0 =
        has_variable_surf_roughness ? &repo.get_field("lowerz0") : nullptr;
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = repo.mesh().Geom(lev);
        const auto& domain = geom.Domain();
        // Read dz
        const amrex::Real dz = geom.CellSize()[idim];
        amrex::MFItInfo mfi_info{};
        const auto& rho_lev = density(lev);
        const auto& vold_lev = velocity.state(FieldState::Old)(lev);
        const auto& told_lev = temperature.state(FieldState::Old)(lev);
        auto& theta = temperature(lev);
        const auto& eta_lev = alpha(lev);

        if (amrex::Gpu::notInLaunchRegion()) {
            mfi_info.SetDynamic(true);
        }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(theta, mfi_info); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.validbox();
            const auto& vold_arr = vold_lev.const_array(mfi);
            const auto& told_arr = told_lev.const_array(mfi);
            const auto& tarr = theta.array(mfi);
            const auto& den = rho_lev.const_array(mfi);
            const auto& eta = eta_lev.const_array(mfi);
            const auto& z0_arr = (has_variable_surf_roughness)
                                     ? (*m_terrainz0)(lev).const_array(mfi)
                                     : amrex::Array4<double>();
            const auto& surf_temp_arr =
                (has_variable_surf_temp) ? (*m_surf_temp)(lev).const_array(mfi)
                                         : amrex::Array4<double>();
            const amrex::Real louis_bh = m_louis_bh;
            const amrex::Real louis_dh = m_louis_dh;
            if (bx.smallEnd(idim) == domain.smallEnd(idim) &&
                temperature.bc_type()[zlo] == BC::wall_model) {
                amrex::ParallelFor(
                    amrex::bdryLo(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real alphaT = eta(i, j, k);
                        const amrex::Real uu = vold_arr(i, j, k, 0);
                        const amrex::Real vv = vold_arr(i, j, k, 1);
                        const amrex::Real wspd = std::sqrt(uu * uu + vv * vv);
                        const amrex::Real theta2 = told_arr(i, j, k);
                        amrex::Real ustarthetastar = 0.0;
                        if (has_variable_surf_roughness &&
                            has_variable_surf_temp) {
                            const amrex::Real theta_surface =
                                surf_temp_arr(i, j, k, 1);
                            const amrex::Real rib =
                                (9.81 * dz * (theta2 - theta_surface)) /
                                (theta2 * wspd * wspd);
                            const amrex::Real louis_ah =
                                0.41 / std::log(0.5 * dz / z0_arr(i, j, k));
                            const amrex::Real louis_ch =
                                5.3 * std::pow(louis_ah, 2) * louis_bh *
                                std::sqrt(0.5 * dz / z0_arr(i, j, k));
                            const amrex::Real Fh =
                                (rib > 0)
                                    ? 1 / std::pow(1 + louis_dh * rib, 2)
                                    : 1 - louis_bh * rib /
                                              (1 +
                                               louis_ch *
                                                   std::sqrt(std::abs(rib)));
                            ustarthetastar = 1.0 / 0.7 * std::pow(louis_ah, 2) *
                                             wspd * (theta2 - theta_surface) *
                                             Fh;
                        }
                        tarr(i, j, k - 1) =
                            (has_variable_surf_roughness)
                                ? den(i, j, k) / alphaT * ustarthetastar
                                : den(i, j, k) * tau.calc_theta(wspd, theta2) /
                                      alphaT;
                    });
            }
        }
    }
}

void ABLTempWallFunc::operator()(Field& temperature, const FieldState rho_state)
{

    const auto& mo = m_wall_func.mo();

    if (m_wall_shear_stress_type == "moeng") {

        auto tau = ShearStressMoeng(mo);
        wall_model(temperature, rho_state, tau);

    } else if (m_wall_shear_stress_type == "constant") {

        auto tau = ShearStressConstant(mo);
        wall_model(temperature, rho_state, tau);

    } else if (m_wall_shear_stress_type == "local") {

        auto tau = ShearStressLocal(mo);
        wall_model(temperature, rho_state, tau);

    } else if (m_wall_shear_stress_type == "schumann") {

        auto tau = ShearStressSchumann(mo);
        wall_model(temperature, rho_state, tau);

    } else if (m_wall_shear_stress_type == "donelan") {

        auto tau = ShearStressDonelan(mo);
        wall_model(temperature, rho_state, tau);
    }
}

} // namespace amr_wind
