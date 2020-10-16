#ifndef MODATA_H
#define MODATA_H

#include "amr-wind/utilities/trig_ops.H"

namespace amr_wind {

/** Monin-Obukhov surface layer profile
 *  \ingroup we_abl
 *
 * van der Laan, P., Kelly, M. C., & Sørensen, N. N. (2017). A new k-epsilon
 * model consistent with Monin-Obukhov similarity theory. Wind Energy,
 * 20(3), 479–489. https://doi.org/10.1002/we.2017
 *
 * Consistent with Dyer (1974) formulation from page 57, Chapter 2, Modeling
 * the vertical ABL structure in Modelling of Atmospheric Flow Fields,
 * Demetri P Lalas and Corrado F Ratto, January 1996,
 * https://doi.org/10.1142/2975.
 */
struct MOData
{
    enum ThetaCalcType {
        HEAT_FLUX = 0,      ///< Heat-flux specified
        SURFACE_TEMPERATURE ///< Surface temperature specified
    };

    amrex::Real zref{0.0};     ///< Reference height (m)
    amrex::Real z0{0.1};       ///< Roughness height (m)
    amrex::Real utau;          ///< Friction velocity (m/s)
    amrex::Real kappa{0.41};   ///< von Karman constant
    amrex::Real gravity{9.81}; ///< Acceleration due to gravity (m/s^2)
    amrex::Real obukhov_len{1.0e16}; ///< Non-dimensional Obukhov length

    amrex::Real vel_mean[AMREX_SPACEDIM]; ///< Mean velocity (at zref)
    amrex::Real vmag_mean;                ///< Mean wind speed (at zref)
    amrex::Real theta_mean;               ///< Mean potential temperature

    amrex::Real surf_temp_flux; ///< Heat flux
    amrex::Real surf_temp;      ///< Instantaneous surface temperature
    amrex::Real ref_temp;       ///< Reference temperature

    amrex::Real gamma_m{5.0};
    amrex::Real gamma_h{5.0};
    amrex::Real beta_m{16.0};

    ThetaCalcType alg_type{HEAT_FLUX};

    amrex::Real phi_m() const
    {
        return std::log(zref / z0) - calc_psi_m(zref / obukhov_len);
    }

    amrex::Real phi_m(amrex::Real z) const
    {
        return std::log(z / z0) - calc_psi_m(z / obukhov_len);
    }

    amrex::Real phi_h() const
    {
        return std::log(zref / z0) - calc_psi_h(zref / obukhov_len);
    }

    amrex::Real phi_h(amrex::Real z) const
    {
        return std::log(z / z0) - calc_psi_h(z / obukhov_len);
    }

    amrex::Real calc_psi_m(amrex::Real zeta) const;
    amrex::Real calc_psi_h(amrex::Real zeta) const;
    void update_fluxes(int max_iters = 25);
};

} // namespace amr_wind

#endif /* MODATA_H */
