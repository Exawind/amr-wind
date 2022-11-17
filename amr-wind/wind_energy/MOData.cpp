#include <limits>
#include "amr-wind/wind_energy/MOData.H"

namespace amr_wind {

/*
 * van der Laan, P., Kelly, M. C., & Sørensen, N. N. (2017). A new k-epsilon
 * model consistent with Monin-Obukhov similarity theory. Wind Energy,
 * 20(3), 479–489. https://doi.org/10.1002/we.2017
 *
 * Consistent with Dyer (1974) formulation from page 57, Chapter 2, Modeling
 * the vertical ABL structure in Modelling of Atmospheric Flow Fields,
 * Demetri P Lalas and Corrado F Ratto, January 1996,
 * https://doi.org/10.1142/2975.
 */

amrex::Real MOData::calc_psi_m(amrex::Real zeta) const
{
    if (zeta > 0) {
        return -gamma_m * zeta;
    }
    amrex::Real x = std::sqrt(std::sqrt(1 - beta_m * zeta));
    return 2.0 * std::log(0.5 * (1.0 + x)) + log(0.5 * (1 + x * x)) -
           2.0 * std::atan(x) + utils::half_pi();
}

amrex::Real MOData::calc_psi_h(amrex::Real zeta) const
{
    if (zeta > 0) {
        return -gamma_h * zeta;
    }
    amrex::Real x = std::sqrt(1 - beta_h * zeta);
    return 2.0 * std::log(0.5 * (1 + x));
}

void MOData::update_fluxes(int max_iters)
{
    constexpr amrex::Real eps = 1.0e-16;
    amrex::Real zeta = 0.0;
    amrex::Real utau_iter = 0.0;

    // Initialize variables
    amrex::Real psi_m = 0.0;
    amrex::Real psi_h = 0.0;
    utau = kappa * vmag_mean / (std::log(zref / z0));

    std::cout << "update_fluxes():" << std::endl;
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

        if (std::abs(surf_temp_flux) > eps) {
            // Stable and unstable ABL conditions
            obukhov_len = -utau * utau * utau * theta_mean /
                          (kappa * gravity * surf_temp_flux);
            zeta = zref / obukhov_len;
        } else {
            // Neutral conditions
            obukhov_len = std::numeric_limits<amrex::Real>::max();
            zeta = 0.0;
        }
        psi_m = calc_psi_m(zeta);
        psi_h = calc_psi_h(zeta);
        utau = kappa * vmag_mean / (std::log(zref / z0) - psi_m);
/*
        std::cout << "    -iter: " << iter << std::endl;
        std::cout << "        utau = " << utau << std::endl;
        std::cout << "        kappa = " << kappa << std::endl;
        std::cout << "        vmag_mean = " << vmag_mean << std::endl;
        std::cout << "        zref = " << zref << std::endl;
        std::cout << "        z0 = " << z0 << std::endl;
        std::cout << "        psi_m = " << psi_m << std::endl;
        ++iter;
*/
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

} // namespace amr_wind
