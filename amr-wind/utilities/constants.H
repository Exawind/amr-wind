#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "AMReX_REAL.H"
#include "AMReX_Math.H"
#include "AMReX_GpuQualifiers.H"

namespace amr_wind::constants {

//! A number close to zero
static constexpr amrex::Real SMALL_NUM = static_cast<amrex::Real>(
    std::numeric_limits<amrex::Real>::epsilon() * 1e10);

//! A number very close to zero
static constexpr amrex::Real EPS =
    static_cast<amrex::Real>(std::numeric_limits<amrex::Real>::epsilon());

//! A tight tolerance
static constexpr amrex::Real TIGHT_TOL = 1e-12;

//! A loose tolerance
static constexpr amrex::Real LOOSE_TOL = 1e-8;

//! A large negative number
static constexpr amrex::Real LOW_NUM = static_cast<amrex::Real>(
    std::numeric_limits<amrex::Real>::lowest() * 1e-10);

//! A large positive number
static constexpr amrex::Real LARGE_NUM =
    static_cast<amrex::Real>(std::numeric_limits<amrex::Real>::max() * 1e-10);

//! Closeness comparison
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE bool
is_close(const amrex::Real a, const amrex::Real b)
{
    return amrex::Math::abs(a - b) < EPS;
}

} // namespace amr_wind::constants
#endif
