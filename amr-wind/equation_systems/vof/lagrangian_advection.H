#ifndef LAGRANGIAN_ADVECTION_H_
#define LAGRANGIAN_ADVECTION_H_

#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include <AMReX_FArrayBox.H>
#include <cmath>

namespace amr_wind {
namespace multiphase {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void lagrangian_advection(
    const int i,
    const int j,
    const int k,
    const int dim,
    const amrex::Real dtdx,
    const amrex::Real velL,
    const amrex::Real velR,
    amrex::Array4<amrex::Real> const& volfrac,
    amrex::Array4<amrex::Real> const& fluxL,
    amrex::Array4<amrex::Real> const& fluxC,
    amrex::Array4<amrex::Real> const& fluxR)
{
    using namespace amrex;

    constexpr Real tiny = 1e-12;
    Real mx = 0.0, my = 0.0, mz = 0.0, alpha = 0.0;
    Real mmL, mmR;
    Real aL = velL * dtdx;
    Real aR = velR * dtdx;

    fluxL(i, j, k) = 0.0;
    fluxC(i, j, k) = 0.0;
    fluxR(i, j, k) = 0.0;

    if (std::abs(volfrac(i, j, k) - 1.0) <= tiny) {
        fluxL(i, j, k) = std::max(-aL, 0.0);
        fluxC(i, j, k) = 1.0 - std::max(aL, 0.0) + std::min(aR, 0.0);
        fluxR(i, j, k) = std::max(aR, 0.0);
    } else if (volfrac(i, j, k) > 0.0) {
        fit_plane(i, j, k, volfrac, mx, my, mz, alpha);
        mmL = std::max(aL, 0.0);
        mmR = 1.0 - mmL + std::min(0.0, aR);
        if (dim == 1) {
            mx = mx / (1.0 - aL + aR);
            alpha = alpha + mx * aL;
            if (aL < 0.0)
                fluxL(i, j, k) = cut_volume(mx, my, mz, alpha, aL, -aL);
            if (aR > 0.0)
                fluxR(i, j, k) = cut_volume(mx, my, mz, alpha, 1.0, aR);
            fluxC(i, j, k) = cut_volume(mx, my, mz, alpha, mmL, mmR);
        } else if (dim == 2) {
            my = my / (1.0 - aL + aR);
            alpha = alpha + my * aL;
            if (aL < 0.0)
                fluxL(i, j, k) = cut_volume(my, mz, mx, alpha, aL, -aL);
            if (aR > 0.0)
                fluxR(i, j, k) = cut_volume(my, mz, mx, alpha, 1.0, aR);
            fluxC(i, j, k) = cut_volume(my, mz, mx, alpha, mmL, mmR);
        } else if (dim == 3) {
            mz = mz / (1.0 - aL + aR);
            alpha = alpha + mz * aL;
            if (aL < 0.0)
                fluxL(i, j, k) = cut_volume(mz, mx, my, alpha, aL, -aL);
            if (aR > 0.0)
                fluxR(i, j, k) = cut_volume(mz, mx, my, alpha, 1.0, aR);
            fluxC(i, j, k) = cut_volume(mz, mx, my, alpha, mmL, mmR);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void balance_fluxes(
    const int i,
    const int j,
    const int k,
    const int dim,
    amrex::Array4<amrex::Real> const& volfrac,
    amrex::Array4<amrex::Real> const& fluxL,
    amrex::Array4<amrex::Real> const& fluxC,
    amrex::Array4<amrex::Real> const& fluxR,
    amrex::BCRec const* bc,
    const int domlo,
    const int domhi)
{

    // Before balance we need to take care of the BC fluxes
    // TO DO

    if (dim == 1) {
        volfrac(i, j, k) =
            fluxR(i - 1, j, k) + fluxC(i, j, k) + fluxL(i + 1, j, k);
    } else if (dim == 2) {
        volfrac(i, j, k) =
            fluxR(i, j - 1, k) + fluxC(i, j, k) + fluxL(i, j + 1, k);
    } else if (dim == 3) {
        volfrac(i, j, k) =
            fluxR(i, j, k - 1) + fluxC(i, j, k) + fluxL(i, j, k + 1);
    }
    // Do clipping
    // volfrac(i, j, k) = std::max(0.0, std::min(1.0, volfrac(i, j, k)));
}
} // namespace multiphase
} // namespace amr_wind
#endif // LAGRANGIAN_ADVECTION.H