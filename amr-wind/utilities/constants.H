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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
is_close(const amrex::Real a, const amrex::Real b)
{
    return amrex::Math::abs(a - b) < EPS;
}

//! Boltzmann constant (J/K)
static constexpr amrex::Real BOLTZMANN_CONSTANT = 1.380649 * 1e-23;

//! Avogadro's constant (mol^-1)
static constexpr amrex::Real AVOGADRO_CONSTANT = 6.02214076 * 1e23;

//! Universal gas constant (J/K mol)
static constexpr amrex::Real UNIVERSAL_GAS_CONSTANT =
    AVOGADRO_CONSTANT * BOLTZMANN_CONSTANT;

//! Heat capacity ratio
static constexpr amrex::Real HEAT_CAPACITY_RATIO = 1.4;

//! Molar mass of air (kg/mol)
static constexpr amrex::Real MOLAR_MASS_AIR = 0.02896492;

//! Molar mass of water vapor (kg/mol)
static constexpr amrex::Real MOLAR_MASS_WATER_VAPOR = 0.01801528;

} // namespace amr_wind::constants
#endif
