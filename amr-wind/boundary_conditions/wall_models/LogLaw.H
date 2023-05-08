#ifndef LogLaw_H
#define LogLaw_H

#include "AMReX_AmrCore.H"
namespace amr_wind {
struct LogLaw
{
    /*
     * A simple wall model that sets the wall-shear stress
     * based on computing u_tau given the horizontal velocity
     * magnitude at a zref. This is akin to an explicit non-linear
     * Robin boundary condition at the wall.
     */

    // Log law constants from Lee & Moser 2015
    // https://doi.org/10.1017/jfm.2015.268.
    amrex::Real B{4.27};
    amrex::Real kappa{0.384};
    int max_iters = 25; // Max iterations for u_tau Newton-Raphson solve
    // Reference height for log law
    amrex::Real zref;
    amrex::Real nu; // molecular viscosity
    // u_tau state variable, gets updated in update_utau
    amrex::Real utau_mean;
    amrex::Real wspd_mean; // mean horizontal velocity magnitude

    void update_utau_mean();
    amrex::Real get_utau(amrex::Real wspd) const;
};
} /* namespace amr_wind */

#endif /* LogLaw_H */
