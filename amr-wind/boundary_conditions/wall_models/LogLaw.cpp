#include "amr-wind/boundary_conditions/wall_models/LogLaw.H"

namespace amr_wind{
amrex::Real LogLaw::get_utau(amrex::Real wspd) const
{
    amrex::Real utau_iter;
    amrex::Real wspd_pred;
    amrex::Real wspd_deriv;
    amrex::Real zplus;
    amrex::Real utau = utau_mean;
    int iter = 0;
    do{
        utau_iter = utau;
        zplus = zref * utau/nu;
        wspd_pred = utau * (std::log(zplus)/kappa + B); //Get wspd for a given utau from log-law
        wspd_deriv = (1 + std::log(zplus))/kappa + B; // d(wspd)/d(utau)
        utau = utau - (wspd_pred - wspd)/wspd_deriv; // Newton-Raphson update
        ++iter;
    } while ((std::abs(utau_iter - utau) > 1e-5) && iter <= max_iters);

    if (iter >= max_iters) {
        amrex::Print()
            << "LogLaw::update_utau: Convergence criteria not met after "
            << max_iters << " iterations" << std::endl
            << "utau = " << utau << " wspd = " << wspd << std::endl;
    }
    return utau;
}
void LogLaw::update_utau_mean()
{
    utau_mean = get_utau(wspd_mean);
}
} // namespace amr_wind
