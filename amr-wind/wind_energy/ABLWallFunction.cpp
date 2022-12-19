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
    amrex::ParmParse pp("ABL");

    pp.query("kappa", m_mo.kappa);
    pp.query("mo_gamma_m", m_mo.gamma_m);
    pp.query("mo_gamma_h", m_mo.gamma_h);
    pp.query("mo_beta_m", m_mo.beta_m);
    pp.query("mo_beta_h", m_mo.beta_h);
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

    if (pp.contains("surface_temp_flux")) {
        pp.query("surface_temp_flux", m_mo.surf_temp_flux);
        if (pp.contains("surface_temp_init")) {
            pp.get("surface_temp_init", m_surf_temp_init);
        } else {
            amrex::Print()
                << "ABLWallFunction: Initial surface temperature not found for "
                   "ABL. Assuming to be equal to the reference temperature "
                << m_mo.ref_temp << std::endl;
            m_surf_temp_init = m_mo.ref_temp;
        }
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

    m_mo.alg_type =
        m_tempflux ? MOData::HEAT_FLUX : MOData::SURFACE_TEMPERATURE;
    m_mo.g = utils::vec_mag(m_gravity.data());

    if (m_tempflux) {
        m_mo.surf_temp = m_surf_temp_init;
    }
}

void ABLWallFunction::init_log_law_height()
{
    if (m_use_fch) {
        const auto& geom = m_mesh.Geom(0);
        m_mo.zref =
            (geom.ProbLo(m_direction) + 0.5 * geom.CellSize(m_direction));
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
                amrex::max(time.current_time() - m_surf_temp_rate_tstart, 0.0) /
                3600.0;
    }

    if (m_inflow_outflow) {
        m_mo.vel_mean[0] = m_wf_vel[0];
        m_mo.vel_mean[1] = m_wf_vel[1];
        m_mo.vmag_mean = m_wf_vmag;
        m_mo.Su_mean = 0.0; // FIXME: need to fill this correctly
        m_mo.Sv_mean = 0.0; // FIXME: need to fill this correctly
        m_mo.Stheta_mean = 0.0; // FIXME: need to fill this correctly
        m_mo.theta_mean = m_wf_theta;
    } else {
        m_mo.vel_mean[0] = vpa.line_average_interpolated(m_mo.zref, 0);
        m_mo.vel_mean[1] = vpa.line_average_interpolated(m_mo.zref, 1);
        m_mo.vmag_mean = vpa.line_hvelmag_average_interpolated(m_mo.zref);
        m_mo.Su_mean = vpa.line_Su_average_interpolated(m_mo.zref);
        m_mo.Sv_mean = vpa.line_Sv_average_interpolated(m_mo.zref);
        m_mo.Stheta_mean = vpa.line_Stheta_average_interpolated(m_mo.zref);
        m_mo.theta_mean = tpa.line_average_interpolated(m_mo.zref, 0);
    }

    m_mo.update_fluxes();
}

ABLVelWallFunc::ABLVelWallFunc(
    Field& /*unused*/, const ABLWallFunction& wall_func)
    : m_wall_func(wall_func)
{
    amrex::ParmParse pp("ABL");
    pp.query("wall_shear_stress_type", m_wall_shear_stress_type);
    m_wall_shear_stress_type = amrex::toLower(m_wall_shear_stress_type);

    if (m_wall_shear_stress_type == "constant" ||
        m_wall_shear_stress_type == "default" ||
        m_wall_shear_stress_type == "local" ||
        m_wall_shear_stress_type == "schumann" ||
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
    const auto& density = repo.get_field("density", rho_state);
    const auto& viscosity = repo.get_field("velocity_mueff");
    const int nlevels = repo.num_active_levels();

    amrex::Orientation zlo(amrex::Direction::z, amrex::Orientation::low);
    amrex::Orientation zhi(amrex::Direction::z, amrex::Orientation::high);
    if (!(velocity.bc_type()[zlo] == BC::wall_model ||
          velocity.bc_type()[zhi] == BC::wall_model)) {
        return;
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = repo.mesh().Geom(lev);
        const auto& domain = geom.Domain();
        amrex::MFItInfo mfi_info{};

        const auto& rho_lev = density(lev);
        auto& vold_lev = velocity.state(FieldState::Old)(lev);
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
            auto varr = vel_lev.array(mfi);
            auto vold_arr = vold_lev.array(mfi);
            auto den = rho_lev.array(mfi);
            auto eta = eta_lev.array(mfi);

            if (bx.smallEnd(idim) == domain.smallEnd(idim) &&
                velocity.bc_type()[zlo] == BC::wall_model) {
                amrex::ParallelFor(
                    amrex::bdryLo(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real mu = eta(i, j, k);
                        const amrex::Real uu = vold_arr(i, j, k, 0);
                        const amrex::Real vv = vold_arr(i, j, k, 1);
                        const amrex::Real wspd = std::sqrt(uu * uu + vv * vv);

                        // Dirichlet BC
                        varr(i, j, k - 1, 2) = 0.0;

                        // Shear stress BC
                        varr(i, j, k - 1, 0) =
                            tau.calc_vel_x(uu, wspd) * den(i, j, k) / mu;
                        varr(i, j, k - 1, 1) =
                            tau.calc_vel_y(vv, wspd) * den(i, j, k) / mu;
                    });
            }

            if (bx.bigEnd(idim) == domain.bigEnd(idim) &&
                velocity.bc_type()[zhi] == BC::wall_model) {
                amrex::ParallelFor(
                    amrex::bdryHi(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real mu = eta(i, j, k - 1);
                        const amrex::Real uu = vold_arr(i, j, k - 1, 0);
                        const amrex::Real vv = vold_arr(i, j, k - 1, 1);
                        const amrex::Real wspd = std::sqrt(uu * uu + vv * vv);

                        // Dirichlet BC
                        varr(i, j, k, 2) = 0.0;

                        // Shear stress BC
                        varr(i, j, k, 0) =
                            -tau.calc_vel_x(uu, wspd) * den(i, j, k - 1) / mu;
                        varr(i, j, k, 1) =
                            -tau.calc_vel_y(vv, wspd) * den(i, j, k - 1) / mu;
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
    } else {

        auto tau = ShearStressDefault(mo);
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
    amrex::Print() << "Heat Flux model: " << m_wall_shear_stress_type
                   << std::endl;
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

    if (!(temperature.bc_type()[zlo] == BC::wall_model ||
          temperature.bc_type()[zhi] == BC::wall_model)) {
        return;
    }

    BL_PROFILE("amr-wind::ABLTempWallFunc");
    auto& velocity = repo.get_field("velocity");
    const auto& density = repo.get_field("density", rho_state);
    const auto& alpha = repo.get_field("temperature_mueff");
    const int nlevels = repo.num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = repo.mesh().Geom(lev);
        const auto& domain = geom.Domain();
        amrex::MFItInfo mfi_info{};

        const auto& rho_lev = density(lev);
        auto& vold_lev = velocity.state(FieldState::Old)(lev);
        auto& told_lev = temperature.state(FieldState::Old)(lev);
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
            auto vold_arr = vold_lev.array(mfi);
            auto told_arr = told_lev.array(mfi);
            auto tarr = theta.array(mfi);
            auto den = rho_lev.array(mfi);
            auto eta = eta_lev.array(mfi);

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
                        tarr(i, j, k - 1) = den(i, j, k) *
                                            tau.calc_theta(wspd, theta2) /
                                            alphaT;
                    });
            }

            if (bx.bigEnd(idim) == domain.bigEnd(idim) &&
                temperature.bc_type()[zhi] == BC::wall_model) {

                amrex::ParallelFor(
                    amrex::bdryHi(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real alphaT = eta(i, j, k - 1);
                        const amrex::Real uu = vold_arr(i, j, k - 1, 0);
                        const amrex::Real vv = vold_arr(i, j, k - 1, 1);
                        const amrex::Real wspd = std::sqrt(uu * uu + vv * vv);
                        const amrex::Real theta2 = told_arr(i, j, k - 1);
                        tarr(i, j, k) = -den(i, j, k - 1) *
                                        tau.calc_theta(wspd, theta2) / alphaT;
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
    } else {

        auto tau = ShearStressDefault(mo);
        wall_model(temperature, rho_state, tau);
    }
}

} // namespace amr_wind
