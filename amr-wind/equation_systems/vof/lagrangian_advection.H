#ifndef LAGRANGIAN_ADVECTION_H_
#define LAGRANGIAN_ADVECTION_H_

#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
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
    int bclo = bc->lo(dim - 1);
    int bchi = bc->hi(dim - 1);
    /*
        amrex::Print() << dim << " " << bclo << " " << bchi << std::endl;
        amrex::Print() << domlo << " " << domhi << std::endl;
    */
    if (dim == 1) {
        if (bclo == amrex::BCType::ext_dir || bclo == amrex::BCType::hoextrap) {
            if (i == domlo) {
                fluxR(domlo - 1, j, k) = 0.0;
            }
        }
        if (bchi == amrex::BCType::ext_dir || bchi == amrex::BCType::hoextrap) {
            if (i == domhi) {
                fluxL(domhi + 1, j, k) = 0.0;
            }
        }
        volfrac(i, j, k) =
            fluxR(i - 1, j, k) + fluxC(i, j, k) + fluxL(i + 1, j, k);
    } else if (dim == 2) {
        if (bclo == amrex::BCType::ext_dir || bclo == amrex::BCType::hoextrap) {
            if (j == domlo) {
                fluxR(i, domlo - 1, k) = 0.0;
            }
        }
        if (bchi == amrex::BCType::ext_dir || bchi == amrex::BCType::hoextrap) {
            if (j == domhi) {
                fluxL(i, domhi + 1, k) = 0.0;
            }
        }
        volfrac(i, j, k) =
            fluxR(i, j - 1, k) + fluxC(i, j, k) + fluxL(i, j + 1, k);
    } else if (dim == 3) {
        if (bclo == amrex::BCType::ext_dir || bclo == amrex::BCType::hoextrap) {
            if (k == domlo) {
                fluxR(i, j, domlo - 1) = 0.0;
            }
        }
        if (bchi == amrex::BCType::ext_dir || bchi == amrex::BCType::hoextrap) {
            if (k == domhi) {
                fluxL(i, j, domhi + 1) = 0.0;
            }
        }
        volfrac(i, j, k) =
            fluxR(i, j, k - 1) + fluxC(i, j, k) + fluxL(i, j, k + 1);
    }
    // Do clipping
    volfrac(i, j, k) = std::max(0.0, std::min(1.0, volfrac(i, j, k)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void remove_vof_debris(
    const int i,
    const int j,
    const int k,
    amrex::Array4<amrex::Real> const& volfrac)
{
    amrex::Real small_vof = 1e-6;
    amrex::Real volFxL = volfrac(i - 1, j, k);
    amrex::Real volFxR = volfrac(i + 1, j, k);
    amrex::Real volFyL = volfrac(i, j - 1, k);
    amrex::Real volFyR = volfrac(i, j + 1, k);
    amrex::Real volFzL = volfrac(i, j, k - 1);
    amrex::Real volFzR = volfrac(i, j, k + 1);

    if (volfrac(i, j, k) > 0.0 && volFxL < small_vof && volFxR < small_vof &&
        volFyL < small_vof && volFyR < small_vof && volFzL < small_vof &&
        volFzR < small_vof) {
        volfrac(i, j, k) = 0.0;
    }
}

} // namespace multiphase
} // namespace amr_wind
#endif // LAGRANGIAN_ADVECTION.H