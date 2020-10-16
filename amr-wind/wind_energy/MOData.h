#ifndef MODATA_H
#define MODATA_H

#include "amr-wind/utilities/trig_ops.H"

namespace amr_wind {

struct MOData
{
    enum ThetaCalcType {
        HEAT_FLUX = 0,
        SURFACE_TEMPERATURE
    };

    amrex::Real zref{0.0};     ///< Reference height (m)
    amrex::Real z0{0.1};       ///< Roughness height (m)
    amrex::Real utau;          ///< Friction velocity (m/s)
    amrex::Real kappa{0.41};   ///< von Karman constant
    amrex::Real gravity{9.81}; ///< Acceleration due to gravity (m/s^2)
    amrex::Real obukhov_len{1.0e16};

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

    amrex::Real calc_psi_m(amrex::Real zeta);
    amrex::Real calc_psi_h(amrex::Real zeta);
    void update_fluxes(int max_iters = 25);
};

} // namespace amr_wind

#endif /* MODATA_H */
