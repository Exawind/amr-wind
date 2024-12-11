#ifndef TRIG_OPS_H
#define TRIG_OPS_H

/** \file trig_ops.H
 *
 *  Trigonometric utilities
 */

#include <cmath>
#include "AMReX_REAL.H"
#include "AMReX_Gpu.H"

namespace amr_wind::utils {

//! Return \f$\pi\f$ as an amrex::Real
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr amrex::Real pi()
{
    return M_PI;
}

//! Return \f$2 \pi\f$
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr amrex::Real two_pi()
{
    return 2.0 * M_PI;
}

//! Return \f$\pi / 2\f$
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr amrex::Real half_pi()
{
    return 0.5 * M_PI;
}

//! Convert from degrees to radians
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
radians(const amrex::Real deg_val)
{
    return pi() * deg_val / 180.0;
}

//! Convert from radians to degrees
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
degrees(const amrex::Real rad_val)
{
    return 180.0 * rad_val / pi();
}

} // namespace amr_wind::utils

#endif /* TRIG_OPS_H */
