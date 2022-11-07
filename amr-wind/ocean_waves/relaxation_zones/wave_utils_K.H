#ifndef WAVES_UTILS_H_
#define WAVES_UTILS_H_

#include <AMReX_FArrayBox.H>
#include <cmath>

namespace amr_wind {
namespace ocean_waves {
namespace relaxation_zones {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
free_surface_to_vof(amrex::Real eta, amrex::Real z, amrex::Real dz)
{
    amrex::Real volfrac = 0.0;
    if (eta - z >= dz / 2.) {
        volfrac = 1.;
    } else if (eta - z <= -dz / 2.) {
        volfrac = 0.;
    } else if (std::abs(eta - z) < dz / 2.) {
        if (eta <= z) {
            volfrac = 1. - (z - eta + dz / 2.) / dz;
        } else {
            volfrac = (eta - z + dz / 2.) / dz;
        }
    }

    return volfrac;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
Gamma_generate(amrex::Real x, amrex::Real gen_length)
{
    amrex::Real Gamma = 0.0;
    amrex::Real xtilde = std::max(std::min(1. - x / gen_length, 1.0), 0.0);
    Gamma = 1.0 - std::expm1(std::pow(xtilde, 3.5)) / std::expm1(1.0);
    return Gamma;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real Gamma_absorb(
    amrex::Real x, amrex::Real absorb_length, amrex::Real absorb_length_factor)
{
    amrex::Real Gamma = 0.0;
    amrex::Real xtilde =
        std::max(std::abs(x) / (absorb_length * absorb_length_factor), 0.0);
    Gamma = 1.0 - std::expm1(std::pow(xtilde, 3.5)) / std::expm1(1.0);
    return Gamma;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
ramp(amrex::Real time, amrex::Real ramp_period)
{
    amrex::Real f = 1.0;
    if (time < ramp_period) {
        f = time / ramp_period -
            (1.0 / M_PI) * std::sin(M_PI * (time / ramp_period));
    }
    return f;
}

} // namespace relaxation_zones
} // namespace ocean_waves
} // namespace amr_wind

#endif
