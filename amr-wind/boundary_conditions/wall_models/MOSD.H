#ifndef MOSD_H
#define MOSD_H

#include "AMReX_AmrCore.H"
#include <AMReX.H>
namespace amr_wind {
struct MOSD
{
    /*
     * A dynamic wall model that calculates the stress from wave to wind
     * based on geometric information of the wave.
     */

    amrex::Real amplitude;
    amrex::Real wavenumber;
    amrex::Real omega;
    amrex::Real time;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    get_dyn_tau(amrex::Real u_dx, amrex::Real v_dx, amrex::Real xc, amrex::Real unit_nor) const
    {
        amrex::Real tau_wave;

        // Building the wave surface, gradients, wave velocities and unit normal
        amrex::Real dx_eta_wave = -amplitude * wavenumber * std::sin(wavenumber * xc  - omega * time);
        amrex::Real dy_eta_wave = 0;
        amrex::Real dt_eta_wave = amplitude * omega * std::sin(wavenumber * xc  - omega * time);
        amrex::Real grad_eta_wave = std::sqrt(dx_eta_wave*dx_eta_wave + dy_eta_wave*dy_eta_wave);
        amrex::Real Cx_wave = - dt_eta_wave * dx_eta_wave *(1/(grad_eta_wave*grad_eta_wave));
        amrex::Real Cy_wave = - dt_eta_wave * dy_eta_wave *(1/(grad_eta_wave*grad_eta_wave));
        amrex::Real n_x =  dx_eta_wave/grad_eta_wave;
        amrex::Real n_y =  dy_eta_wave/grad_eta_wave;
                    
        // Calculating the relative velocity, heaviside function
        amrex::Real u_r = u_dx - Cx_wave;
        amrex::Real v_r = v_dx - Cy_wave;
        amrex::Real ur_mag = std::sqrt(u_r*u_r*n_x*n_x + v_r*v_r*n_y*n_y);
        amrex::Real Heavi_arg = (u_r*dx_eta_wave + v_r*dy_eta_wave);
        amrex::Real Heavi = (Heavi_arg + std::abs(Heavi_arg))/(2*Heavi_arg);

        // Calculating the magnitude of the stress
        if (unit_nor == 0){
        tau_wave =  (1/M_PI)*ur_mag*ur_mag*grad_eta_wave*grad_eta_wave*Heavi*n_x ;
        } else {
        tau_wave =  (1/M_PI)*ur_mag*ur_mag*grad_eta_wave*grad_eta_wave*Heavi*n_y ;
        }
        return tau_wave;
    }
};

} // namespace amr_wind
#endif /* MOSD_H */
