#ifndef WAVES_UTILS_H_
#define WAVES_UTILS_H_

#include <AMReX_FArrayBox.H>
#include <cmath>

using WaveVec = amrex::GpuArray<amrex::Real, 4>;

namespace amr_wind::ocean_waves::utils {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real free_surface_to_vof(
    const amrex::Real eta, const amrex::Real z, const amrex::Real dz)
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
gamma_generate(const amrex::Real x, const amrex::Real gen_length)
{
    const amrex::Real xtilde =
        std::max(std::min(1. - x / gen_length, 1.0), 0.0);
    return (1.0 - std::expm1(std::pow(xtilde, 3.5)) / std::expm1(1.0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real gamma_absorb(
    const amrex::Real x,
    const amrex::Real absorb_length,
    const amrex::Real absorb_length_factor)
{
    const amrex::Real xtilde = std::max(
        std::min(x / (absorb_length * absorb_length_factor), 1.0), 0.0);
    return (1.0 - std::expm1(std::pow(xtilde, 3.5)) / std::expm1(1.0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
ramp(const amrex::Real time, const amrex::Real ramp_period)
{
    amrex::Real f = 1.0;
    if (time < ramp_period) {
        f = time / ramp_period -
            (1.0 / M_PI) * std::sin(M_PI * (time / ramp_period));
    }
    return f;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real combine_linear(
    const amrex::Real factor,
    const amrex::Real target,
    const amrex::Real current)
{
    return (1. - factor) * target + factor * current;
}

} // namespace amr_wind::ocean_waves::utils

#endif
