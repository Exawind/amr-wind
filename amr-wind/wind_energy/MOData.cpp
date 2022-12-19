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


//Non-dimensional velocity shear, phi_m
amrex::Real MOData::phi_m(amrex::Real z) const
{
    // if neutral:
    amrex::Real phi_m = 1.0;

    // if stable:
    if ((L < L_neutral) && (L > 0.0))
    {
        phi_m = 1.0 + gamma_m*(z/L);
    }

    // if unstable:
    else if ((L > -L_neutral) && (L <= 0.0))
    {
        phi_m = std::pow(1.0 - beta_m*(z/L), -0.25);
    }

    return phi_m;
}
amrex::Real MOData::phi_m() const
{
    return phi_m(zref);
}

// Non-dimensional potential temperature shear, phi_h
amrex::Real MOData::phi_h(amrex::Real z) const
{
    // if neutral:
    amrex::Real phi_h = 1.0;

    // if stable:
    if ((L < L_neutral) && (L > 0.0))
    {
        phi_h = alpha_h + gamma_h*(z/L);
    }

    // if unstable:    
    else if ((L > -L_neutral) && (L <= 0.0))
    {
        phi_h = alpha_h*std::pow(1.0 - beta_h*(z/L), -0.5);
    }

    return phi_h;
}
amrex::Real MOData::phi_h() const
{
    return phi_h(zref);
}

// Velocity profile stability correction function, psi_m
amrex::Real MOData::psi_m(amrex::Real z) const
{
    // if neutral:
    amrex::Real psi_m = 0.0;

    // if stable:
    if ((L < L_neutral) && (L > 0.0))
    {
        psi_m = 1.0 - phi_m(z);
    }

    // if unstable:
    else if ((L > -L_neutral) && (L <= 0.0))
    {
        psi_m =   2.0 * std::log((1.0 + std::pow(phi_m(z),-1.0)) / 2.0)
                + 1.0 * std::log((1.0 + std::pow(phi_m(z),-2.0)) / 2.0)
                - 2.0 * std::atan(1.0 / phi_m(z))
                + utils::half_pi();
    }

    return psi_m;    
}

// Potential temperature profile stability correction function, psi_h
amrex::Real MOData::psi_h(amrex::Real z) const
{
    // if neutral:
    amrex::Real psi_h = 0.0;

    // if stable:
    if ((L < L_neutral) && (L > 0.0))
    {
        psi_h = 1.0 - (phi_h(z)/alpha_h);
    }
     //           1.0 - (phi_h_/self.alpha_h)

    // if unstable:
    else if ((L > -L_neutral) && (L <= 0.0))
    {
        psi_h = 2.0 * std::log((1.0 + std::pow(phi_h(z)/alpha_h,-1.0)) / 2.0);
    }

    return psi_h;
}

// Update Obukhov length
void MOData::update_L()
{
    // if non-neutral
    if (std::abs(surf_temp_flux) > 1.0E-6)
    {
        L = (-std::pow(utau,3.0) * surf_temp) / (kappa * g * surf_temp_flux);
    }
    // if neutral
    else
    {
        L = L_neutral;
    }
}

// Update friction velocity
void MOData::update_utau()
{
    utau = (kappa * vmag_mean) / (std::log(zref/z0) - (psi_m(zref) - psi_m(z0)));
    utau = std::max(1.0E-6,utau);
}

// Update the surface temperature
void MOData::update_surf_temp()
{
    surf_temp = theta_mean + 
      ((alpha_h * surf_temp_flux)/(kappa * utau)) * (std::log(zref/z0) - (psi_h(zref) - psi_h(z0)));
}

// Update the surface temperature flux
void MOData::update_surf_temp_flux()
{
    surf_temp_flux = (-(theta_mean - surf_temp) * kappa * utau) /
                           (alpha_h * (std::log(zref/z0) - (psi_h(zref) - psi_h(z0))));
}

void MOData::update_fluxes(int max_iters)
{
    amrex::Real utau1 = std::max(1.0E-6, (kappa * vmag_mean) / (std::log(zref/z0)));
    amrex::Real utau0 = utau1 + 2.0*tol;
    utau = utau1;

    int iter = 0;
    while ((std::abs(utau1 - utau0) > tol) && (iter < max_iters))
    {
        utau0 = utau1;
        update_L();
        update_utau();
        switch (alg_type) 
        {
            case HEAT_FLUX:
                update_surf_temp();
                break;

            case SURFACE_TEMPERATURE:
                update_surf_temp_flux();
                break;
        }
        utau1 = utau;
        iter++;

        std::cout << "    -iter: " << iter << std::endl;
        std::cout << "        utau = " << utau << std::endl;
        std::cout << "        L = " << L << std::endl;
        std::cout << "        kappa = " << kappa << std::endl;
        std::cout << "        vmag_mean = " << vmag_mean << std::endl;
        std::cout << "        zref = " << zref << std::endl;
        std::cout << "        z0 = " << z0 << std::endl;
        std::cout << "        theta_mean = " << theta_mean << std::endl;
        std::cout << "        surf_temp_flux = " << surf_temp_flux << std::endl;
        std::cout << "        surf_temp = " << surf_temp << std::endl;
        std::cout << "        phi_m = " << phi_m() << std::endl;
        std::cout << "        psi_m = " << psi_m(zref) << std::endl;
        std::cout << "        phi_h = " << phi_h() << std::endl;
        std::cout << "        psi_h = " << psi_h(zref) << std::endl;
    }
    



/*
    constexpr amrex::Real eps = 1.0e-16;
    amrex::Real zeta = 0.0;
    amrex::Real utau_iter = 0.0;

    // Initialize variables
    amrex::Real psi_m = 0.0;
    amrex::Real psi_h = 0.0;
    utau = kappa * vmag_mean / (std::log(zref / z0));
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

        std::cout << "    -iter: " << iter << std::endl;
        std::cout << "        utau = " << utau << std::endl;
        std::cout << "        kappa = " << kappa << std::endl;
        std::cout << "        vmag_mean = " << vmag_mean << std::endl;
        std::cout << "        zref = " << zref << std::endl;
        std::cout << "        z0 = " << z0 << std::endl;
        std::cout << "        psi_m = " << psi_m << std::endl;
        ++iter;
    } while ((std::abs(utau_iter - utau) > 1e-5) && iter <= max_iters);
*/

    if (iter >= max_iters) {
        amrex::Print()
            << "MOData::update_fluxes: Convergence criteria not met after "
            << max_iters << " iterations"
            << "\nObuhov length = " << L << " zeta = " << zref/L
            << "\npsi_m = " << psi_m(zref) << " psi_h = " << psi_h(zref)
            << "\nutau = " << utau << " Tsurf = " << surf_temp
            << " q = " << surf_temp_flux << std::endl;
    }
}

} // namespace amr_wind
