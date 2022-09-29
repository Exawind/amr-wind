#ifndef LAGRANGIAN_ADVECTION_H_
#define LAGRANGIAN_ADVECTION_H_

#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
#include <cmath>

namespace amr_wind {
namespace multiphase {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eulerian_implicit(
    const int i,
    const int j,
    const int k,
    const int dir,
    const amrex::Real dtdx,
    const amrex::Real velL,
    const amrex::Real velR,
    amrex::Array4<amrex::Real const> const& volfrac,
    amrex::Array4<amrex::Real> const& fluxL,
    amrex::Array4<amrex::Real> const& fluxR,
    amrex::Array4<amrex::Real> const& vofL,
    amrex::Array4<amrex::Real> const& vofR)
{
    using namespace amrex;

    constexpr Real tiny = 1e-12;
    Real mx = 0.0, my = 0.0, mz = 0.0, alpha = 0.0;
    Real x0, deltax;
    Real aL = velL * dtdx;
    Real aR = velR * dtdx;

    fluxL(i, j, k) = 0.0;
    fluxR(i, j, k) = 0.0;
    vofL(i, j, k) = 0.0;
    vofR(i, j, k) = 0.0;

    if (std::abs(volfrac(i, j, k) - 1.0) <= tiny) {
        fluxL(i, j, k) = amrex::max(-aL, 0.0);
        fluxR(i, j, k) = amrex::max(aR, 0.0);
        if (aL < 0.0) vofL(i, j, k) = 1.0;
        if (aR > 0.0) vofR(i, j, k) = 1.0;
    } else if (volfrac(i, j, k) > tiny) {
        fit_plane(i, j, k, volfrac, mx, my, mz, alpha);
        // Eulerian advection
        x0 = 0.0;
        if (dir == 0) {
            if (aL < 0.0) {
                deltax = -aL;
                vofL(i, j, k) = cut_volume(mx, my, mz, alpha, x0, deltax);
                fluxL(i, j, k) = vofL(i, j, k) * deltax;
            }
            if (aR > 0.0) {
                x0 = 1.0 - aR;
                deltax = aR;
                vofR(i, j, k) = cut_volume(mx, my, mz, alpha, x0, deltax);
                fluxR(i, j, k) = vofR(i, j, k) * deltax;
            }
        } else if (dir == 1) {
            if (aL < 0.0) {
                deltax = -aL;
                vofL(i, j, k) = cut_volume(my, mz, mx, alpha, x0, deltax);
                fluxL(i, j, k) = vofL(i, j, k) * deltax;
            }
            if (aR > 0.0) {
                x0 = 1.0 - aR;
                deltax = aR;
                vofR(i, j, k) = cut_volume(my, mz, mx, alpha, x0, deltax);
                fluxR(i, j, k) = vofR(i, j, k) * deltax;
            }
        } else if (dir == 2) {
            if (aL < 0.0) {
                deltax = -aL;
                vofL(i, j, k) = cut_volume(mz, mx, my, alpha, x0, deltax);
                fluxL(i, j, k) = vofL(i, j, k) * deltax;
            }
            if (aR > 0.0) {
                x0 = 1.0 - aR;
                deltax = aR;
                vofR(i, j, k) = cut_volume(mz, mx, my, alpha, x0, deltax);
                fluxR(i, j, k) = vofR(i, j, k) * deltax;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void balance_eulerian_fluxes(
    const int i,
    const int j,
    const int k,
    const int dir,
    const amrex::Real dxi,
    const amrex::Real velL,
    const amrex::Real velR,
    amrex::Array4<amrex::Real> const& volfrac,
    amrex::Array4<amrex::Real const> const& fluxF,
    amrex::Array4<amrex::Real const> const& fluxC)
{

    if (dir == 0) {
        volfrac(i, j, k) += (fluxF(i, j, k) - fluxF(i + 1, j, k)) * dxi +
                            fluxC(i, j, k) * (velR - velL);
    } else if (dir == 1) {
        volfrac(i, j, k) += (fluxF(i, j, k) - fluxF(i, j + 1, k)) * dxi +
                            fluxC(i, j, k) * (velR - velL);
    } else if (dir == 2) {
        volfrac(i, j, k) += (fluxF(i, j, k) - fluxF(i, j, k + 1)) * dxi +
                            fluxC(i, j, k) * (velR - velL);
    }
    // Do clipping
    volfrac(i, j, k) = amrex::max(0.0, amrex::min(1.0, volfrac(i, j, k)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void fluxes_bc_save(
    const int i,
    const int j,
    const int k,
    const int dir,
    const amrex::Real dx,
    amrex::Box const& bx,
    amrex::Array4<amrex::Real> const& fluxL,
    amrex::Array4<amrex::Real> const& fluxR,
    amrex::Array4<amrex::Real> const& f_f,
    amrex::Array4<amrex::Real> const& vofL,
    amrex::Array4<amrex::Real> const& vofR,
    amrex::Array4<amrex::Real> const& advalpha_f,
    amrex::BCRec const* bc,
    const int domlo,
    const int domhi)
{
    int bclo = bc->lo(dir);
    int bchi = bc->hi(dir);

    if (dir == 0) {
        if (bclo == amrex::BCType::ext_dir || bclo == amrex::BCType::hoextrap) {
            if (i == domlo) {
                fluxR(domlo - 1, j, k) = 0.0;
                vofR(domlo - 1, j, k) = 0.0;
            }
        }
        if (bchi == amrex::BCType::ext_dir || bchi == amrex::BCType::hoextrap) {
            if (i == domhi) {
                fluxL(domhi + 1, j, k) = 0.0;
                vofL(domhi + 1, j, k) = 0.0;
            }
        }
        f_f(i, j, k) = (fluxR(i - 1, j, k) - fluxL(i, j, k)) * dx;
        advalpha_f(i, j, k) = vofR(i - 1, j, k) + vofL(i, j, k);
        if (!bx.contains(i + 1, j, k)) {
            f_f(i + 1, j, k) = (fluxR(i, j, k) - fluxL(i + 1, j, k)) * dx;
            advalpha_f(i + 1, j, k) = vofR(i, j, k) + vofL(i + 1, j, k);
        }
    } else if (dir == 1) {
        if (bclo == amrex::BCType::ext_dir || bclo == amrex::BCType::hoextrap) {
            if (j == domlo) {
                fluxR(i, domlo - 1, k) = 0.0;
                vofR(i, domlo - 1, k) = 0.0;
            }
        }
        if (bchi == amrex::BCType::ext_dir || bchi == amrex::BCType::hoextrap) {
            if (j == domhi) {
                fluxL(i, domhi + 1, k) = 0.0;
                vofL(i, domhi + 1, k) = 0.0;
            }
        }
        f_f(i, j, k) = (fluxR(i, j - 1, k) - fluxL(i, j, k)) * dx;
        advalpha_f(i, j, k) = vofR(i, j - 1, k) + vofL(i, j, k);
        if (!bx.contains(i, j + 1, k)) {
            f_f(i, j + 1, k) = (fluxR(i, j, k) - fluxL(i, j + 1, k)) * dx;
            advalpha_f(i, j + 1, k) = vofR(i, j, k) + vofL(i, j + 1, k);
        }
    } else if (dir == 2) {
        if (bclo == amrex::BCType::ext_dir || bclo == amrex::BCType::hoextrap) {
            if (k == domlo) {
                fluxR(i, j, domlo - 1) = 0.0;
                vofR(i, j, domlo - 1) = 0.0;
            }
        }
        if (bchi == amrex::BCType::ext_dir || bchi == amrex::BCType::hoextrap) {
            if (k == domhi) {
                fluxL(i, j, domhi + 1) = 0.0;
                vofL(i, j, domhi + 1) = 0.0;
            }
        }
        f_f(i, j, k) = (fluxR(i, j, k - 1) - fluxL(i, j, k)) * dx;
        advalpha_f(i, j, k) = vofR(i, j, k - 1) + vofL(i, j, k);
        if (!bx.contains(i, j, k + 1)) {
            f_f(i, j, k + 1) = (fluxR(i, j, k) - fluxL(i, j, k + 1)) * dx;
            advalpha_f(i, j, k + 1) = vofR(i, j, k) + vofL(i, j, k + 1);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void c_mask(
    const int i,
    const int j,
    const int k,
    amrex::Array4<amrex::Real> const& volfrac,
    amrex::Array4<amrex::Real> const& volfrac_masked)
{
    if (volfrac(i, j, k) > 0.5) {
        volfrac_masked(i, j, k) = 1.0;
    } else {
        volfrac_masked(i, j, k) = 0.0;
    }
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
