#include "amr-wind/wind_energy/MOData.h"

namespace amr_wind {

amrex::Real MOData::calc_psi_m(amrex::Real zeta)
{
    if (zeta > 0) {
        return -gamma_m * zeta;
    } else {
        amrex::Real x = std::sqrt(std::sqrt(1 - beta_m * zeta));
        return 2.0 * std::log(0.5 * (1.0 + x)) + log(0.5 * (1 + x * x)) -
            2.0 * std::atan(x) + utils::half_pi();
    }
}

amrex::Real MOData::calc_psi_h(amrex::Real zeta)
{
    if (zeta > 0) {
        return -gamma_h * zeta;
    } else {
        amrex::Real x = std::sqrt(1 - beta_m * zeta);
        return std::log(0.5 * (1 + x));
    }
}

void MOData::update_fluxes(int max_iters)
{
    amrex::Real zeta = 0.0;
    amrex::Real utau_iter = 0.0;

    // Initialize variables
    amrex::Real psi_m = 0.0;
    amrex::Real psi_h = 0.0;
    utau = kappa * vmag_mean / (std::log(zref / z0) - psi_m);

    int iter = 0;
    do {
        utau_iter = utau;
        switch (alg_type) {
        case HEAT_FLUX:
            surf_temp = surf_temp_flux * (std::log(zref / z0) - psi_h) /
                            (utau * kappa) +
                        theta_mean;
            break;

        case SURFACE_TEMPERATURE:
            surf_temp_flux = -(theta_mean - surf_temp) * utau * kappa /
                             (std::log(zref / z0) - psi_h);
            break;
        }

        obukhov_len = -utau * utau * utau * theta_mean /
                      (kappa * gravity * surf_temp_flux);
        zeta = zref / obukhov_len;
        psi_m = calc_psi_m(zeta);
        psi_h = calc_psi_h(zeta);
        utau = kappa * vmag_mean / (std::log(zref / z0) - psi_m);
        ++iter;
    } while ((std::abs(utau_iter - utau) > 1e-5) && iter <= max_iters);

    if (iter >= max_iters) {
        amrex::Print()
            << "MOData::update_fluxes: Convergence criteria not met after "
            << max_iters << " iterations"
            << "\nObuhov length = " << obukhov_len << " zeta = " << zeta
            << "\npsi_m = " << psi_m << " psi_h = " << psi_h
            << "\nutau = " << utau << " Tsurf = " << surf_temp
            << " q = " << surf_temp_flux << std::endl;
    }
}

}
