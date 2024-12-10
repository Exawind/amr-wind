#ifndef WAVES_UTILS_H_
#define WAVES_UTILS_H_

#include <AMReX_FArrayBox.H>
#include <cmath>

namespace amr_wind::ocean_waves::utils {

using WaveVec = amrex::GpuArray<amrex::Real, 4>;

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real free_surface_to_vof(
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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
gamma_generate(const amrex::Real x, const amrex::Real gen_length)
{
    const amrex::Real xtilde =
        std::max(std::min(1. - x / gen_length, 1.0), 0.0);
    return (1.0 - std::expm1(std::pow(xtilde, 3.5)) / std::expm1(1.0));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real gamma_absorb(
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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real combine_linear(
    const amrex::Real factor,
    const amrex::Real target,
    const amrex::Real current)
{
    return (1. - factor) * target + factor * current;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE WaveVec harmonize_profiles_1d(
    const amrex::Real x,
    const amrex::Real left_bdy,
    const amrex::Real left_length,
    const amrex::Real right_bdy,
    const amrex::Real right_length,
    const WaveVec left,
    const WaveVec bulk,
    const WaveVec right)
{
    // Initialize boundary regions entirely with boundary solutions
    if (x <= left_bdy + left_length) {
        return left;
    }
    if (x + right_length >= right_bdy) {
        return right;
    }
    // Combine solutions in central region using relax zone length scales
    const amrex::Real Gamma_left =
        gamma_generate(x - (left_bdy + left_length), 0.5 * left_length);
    const amrex::Real Gamma_right = gamma_absorb(
        x - (right_bdy - 1.5 * right_length), 0.5 * right_length, 1.0);
    const bool right_inactive = Gamma_right > 1.0 - 1e-12;
    WaveVec combo;
    for (int n = 0; n < 4; ++n) {
        combo[n] = right_inactive
                       ? combine_linear(Gamma_left, left[n], bulk[n])
                       : combine_linear(Gamma_right, right[n], bulk[n]);
    }
    return combo;
}

} // namespace amr_wind::ocean_waves::utils

#endif
