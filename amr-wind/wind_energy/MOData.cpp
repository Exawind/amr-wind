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

amrex::Real MOData::calc_psi_m(
    const amrex::Real zeta, const amrex::Real beta_m, const amrex::Real gamma_m)
{
    if (zeta > 0) {
        return -gamma_m * zeta;
    }
    const amrex::Real x = std::sqrt(std::sqrt(1 - beta_m * zeta));
    return 2.0 * std::log(0.5 * (1.0 + x)) + log(0.5 * (1 + x * x)) -
           2.0 * std::atan(x) + utils::half_pi();
}

amrex::Real MOData::calc_psi_m(amrex::Real zeta) const
{
    return calc_psi_m(zeta, beta_m, gamma_m);
}

amrex::Real MOData::calc_psi_h(
    const amrex::Real zeta, const amrex::Real beta_h, const amrex::Real gamma_h)
{
    if (zeta > 0) {
        return -gamma_h * zeta;
    }
    const amrex::Real x = std::sqrt(1 - beta_h * zeta);
    return 2.0 * std::log(0.5 * (1 + x));
}

amrex::Real MOData::calc_psi_h(amrex::Real zeta) const
{
    return calc_psi_h(zeta, beta_h, gamma_h);
}

void MOData::update_fluxes(int max_iters)
{
    constexpr amrex::Real eps = 1.0e-16;
    amrex::Real zeta = 0.0;
    amrex::Real utau_iter = 0.0;

    // Initialize variables
    amrex::Real psi_m_zref = 0.0;
    amrex::Real psi_h_zref = 0.0;
    utau = kappa * vmag_mean / (std::log(zref / z0));

    int iter = 0;
    do {
        utau_iter = utau;
        switch (alg_type) {
        case ThetaCalcType::HEAT_FLUX:
            surf_temp = alpha_h * surf_temp_flux * (std::log(zref / z0t) - psi_h_zref) /
                            (utau * kappa) +
                        theta_mean;
            break;

        case ThetaCalcType::SURFACE_TEMPERATURE:
            surf_temp_flux = -(theta_mean - surf_temp) * utau * kappa /
                             (alpha_h * (std::log(zref / z0t) - psi_h_zref));
            break;

        case ThetaCalcType::NEAR_SURFACE_TEMPERATURE:
            amrex::Real psi_h_zNearSurf = calc_psi_h(zNearSurf / obukhov_len);
            amrex::Real a11 = 1.0;
            amrex::Real a12 = -(alpha_h / (kappa * utau)) * (std::log(zref / z0t) - psi_h_zref);
            amrex::Real a21 = 1.0;
            amrex::Real a22 = -(alpha_h / (kappa * utau)) * (std::log(zNearSurf / z0t) - psi_h_zNearSurf);

            amrex::Real coeff = 1.0 / (a11*a21 - a12*a21);

            surf_temp = coeff * (a22*theta_mean - a12*near_surf_temp);
            surf_temp_flux = coeff * (-a21*theta_mean + a11*near_surf_temp);

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
        psi_m_zref = calc_psi_m(zeta);
        psi_h_zref = calc_psi_h(zeta);
        utau = kappa * vmag_mean / (std::log(zref / z0) - psi_m_zref);
        ++iter;
    } while ((std::abs(utau_iter - utau) > 1e-5) && iter <= max_iters);

    if (iter >= max_iters) {
        amrex::Print()
            << "MOData::update_fluxes: Convergence criteria not met after "
            << max_iters << " iterations\nObuhov length = " << obukhov_len
            << " zeta = (z_ref/L) = " << zeta << "\npsi_m(z_ref/L) = " << psi_m_zref
            << " psi_h(z_ref/L) = " << psi_h_zref << "\nutau = " << utau
            << " Tsurf = " << surf_temp << " q = " << surf_temp_flux
            << std::endl;
    }
}

} // namespace amr_wind
