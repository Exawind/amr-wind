#ifndef GODUNOV_PPM_H
#define GODUNOV_PPM_H

#include <AMReX_Gpu.H>
#include <AMReX_BCRec.H>

/* This header file contains the inlined __host__ __device__ functions required
   for the scalar advection routines for 3D Godunov. It also contains function
   declarations for controlling host functions.  */

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_trans_xbc(
    const int i,
    const int j,
    const int k,
    const int n,
    const amrex::Array4<const amrex::Real>& s,
    amrex::Real& lo,
    amrex::Real& hi,
    amrex::Real& /* uedge */,
    const int bclo,
    const int bchi,
    const int domlo,
    const int domhi)
{
    using namespace amrex;

    // Low X
    if (i <= domlo) {
        if (bclo == BCType::ext_dir) {
            // IAMR does this but it breaks lo/hi symmetry
            // Real st = (uedge <= small_vel) ? hi : s(domlo-1,j,k,n);
            // So here we do something simpler...
            Real st = s(domlo - 1, j, k, n);
            lo = st;
            hi = st;
        }

        else if (
            bclo == BCType::foextrap || bclo == BCType::hoextrap ||
            bclo == BCType::reflect_even) {
            lo = hi;

        } else if (bclo == BCType::reflect_odd) {
            hi = 0.;
            lo = 0.;
        }
    }

    // High X
    if (i > domhi) {
        if (bchi == BCType::ext_dir) {
            // IAMR does this but it breaks lo/hi symmetry
            // Real st = (uedge >= -small_vel)? lo : s(domhi+1,j,k,n);
            // So here we do something simpler...
            Real st = s(domhi + 1, j, k, n);
            lo = st;
            hi = st;
        }

        else if (
            bchi == BCType::foextrap || bchi == BCType::hoextrap ||
            bchi == BCType::reflect_even) {
            hi = lo;

        } else if (bchi == BCType::reflect_odd) {
            hi = 0.;
            lo = 0.;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_trans_ybc(
    const int i,
    const int j,
    const int k,
    const int n,
    const amrex::Array4<const amrex::Real>& s,
    amrex::Real& lo,
    amrex::Real& hi,
    amrex::Real /* vedge */,
    const int bclo,
    const int bchi,
    const int domlo,
    const int domhi)
{
    using namespace amrex;

    // Low Y
    if (j <= domlo) {
        if (bclo == BCType::ext_dir) {
            // IAMR does this but it breaks lo/hi symmetry
            // Real st = (vedge <= small_vel) ? hi : s(i,domlo-1,k,n);
            // So here we do something simpler...
            Real st = s(i, domlo - 1, k, n);
            lo = st;
            hi = st;
        }

        else if (
            bclo == BCType::foextrap || bclo == BCType::hoextrap ||
            bclo == BCType::reflect_even) {
            lo = hi;

        } else if (bclo == BCType::reflect_odd) {
            hi = 0.;
            lo = 0.;
        }
    }

    // High Y
    if (j > domhi) {
        if (bchi == BCType::ext_dir) {
            // IAMR does this but it breaks lo/hi symmetry
            // Real st = (vedge >= -small_vel)? lo : s(i,domhi+1,k,n);
            // So here we do something simpler...
            Real st = s(i, domhi + 1, k, n);
            lo = st;
            hi = st;
        }

        else if (
            bchi == BCType::foextrap || bchi == BCType::hoextrap ||
            bchi == BCType::reflect_even) {
            hi = lo;

        } else if (bchi == BCType::reflect_odd) {
            hi = 0.;
            lo = 0.;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_trans_zbc(
    const int i,
    const int j,
    const int k,
    const int n,
    const amrex::Array4<const amrex::Real>& s,
    amrex::Real& lo,
    amrex::Real& hi,
    amrex::Real /* wedge */,
    const int bclo,
    const int bchi,
    const int domlo,
    const int domhi)
{
    using namespace amrex;

    // Low Z
    if (k <= domlo) {
        if (bclo == BCType::ext_dir) {
            // IAMR does this but it breaks lo/hi symmetry
            // Real st = (wedge <= small_vel) ? hi : s(i,j,domlo-1,n);
            // So here we do something simpler...
            Real st = s(i, j, domlo - 1, n);
            lo = st;
            hi = st;
        }

        else if (
            bclo == BCType::foextrap || bclo == BCType::hoextrap ||
            bclo == BCType::reflect_even) {
            lo = hi;

        } else if (bclo == BCType::reflect_odd) {
            hi = 0.;
            lo = 0.;
        }
    }

    // High Z
    if (k > domhi) {
        if (bchi == BCType::ext_dir) {
            // IAMR does this but it breaks lo/hi symmetry
            // Real st = (wedge >= -small_vel)? lo : s(i,j,domhi+1,n);
            // So here we do something simpler...
            Real st = s(i, j, domhi + 1, n);
            lo = st;
            hi = st;
        } else if (
            bchi == BCType::foextrap || bchi == BCType::hoextrap ||
            bchi == BCType::reflect_even) {
            hi = lo;

        } else if (bchi == BCType::reflect_odd) {
            hi = 0.;
            lo = 0.;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_cc_xbc_lo(
    const int i,
    const int j,
    const int k,
    const int n,
    const amrex::Array4<const amrex::Real>& s,
    amrex::Real& lo,
    amrex::Real& hi,
    const amrex::Array4<const amrex::Real>& uedge,
    const int bclo,
    const int domlo)
{
    using namespace amrex;

    if (i == domlo) {
        if (bclo == BCType::ext_dir && uedge(i, j, k) >= 0.) {
            hi = s(domlo - 1, j, k, n);
            lo = hi;
        } else if (
            bclo == BCType::foextrap || bclo == BCType::hoextrap ||
            bclo == BCType::reflect_even ||
            (bclo == BCType::ext_dir && uedge(i, j, k) < 0)) {
            lo = hi;

        } else if (bclo == BCType::reflect_odd) {
            hi = 0.;
            lo = hi;
        }
    } else {
        return;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_cc_xbc_hi(
    const int i,
    const int j,
    const int k,
    const int n,
    const amrex::Array4<const amrex::Real>& s,
    amrex::Real& lo,
    amrex::Real& hi,
    const amrex::Array4<const amrex::Real>& uedge,
    const int bchi,
    const int domhi)
{
    using namespace amrex;

    if (i == domhi + 1) {
        if (bchi == BCType::ext_dir && uedge(i, j, k) <= 0.) {
            lo = s(domhi + 1, j, k, n);
            hi = lo;
        } else if (
            bchi == BCType::foextrap || bchi == BCType::hoextrap ||
            bchi == BCType::reflect_even ||
            (bchi == BCType::ext_dir && uedge(i, j, k) > 0)) {
            hi = lo;

        } else if (bchi == BCType::reflect_odd) {
            lo = 0.;
            hi = lo;
        }
    } else {
        return;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_cc_ybc_lo(
    const int i,
    const int j,
    const int k,
    const int n,
    const amrex::Array4<const amrex::Real>& s,
    amrex::Real& lo,
    amrex::Real& hi,
    const amrex::Array4<const amrex::Real>& vedge,
    const int bclo,
    const int domlo)
{
    using namespace amrex;

    if (j == domlo) {
        if (bclo == BCType::ext_dir && vedge(i, j, k) >= 0.) {
            hi = s(i, domlo - 1, k, n);
            lo = hi;
        } else if (
            bclo == BCType::foextrap || bclo == BCType::hoextrap ||
            bclo == BCType::reflect_even ||
            (bclo == BCType::ext_dir && vedge(i, j, k) < 0)) {
            lo = hi;

        } else if (bclo == BCType::reflect_odd) {
            hi = 0.;
            lo = hi;
        }
    } else {
        return;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_cc_ybc_hi(
    const int i,
    const int j,
    const int k,
    const int n,
    const amrex::Array4<const amrex::Real>& s,
    amrex::Real& lo,
    amrex::Real& hi,
    const amrex::Array4<const amrex::Real>& vedge,
    const int bchi,
    const int domhi)
{
    using namespace amrex;

    if (j == domhi + 1) {
        if (bchi == BCType::ext_dir && vedge(i, j, k) <= 0.) {
            lo = s(i, domhi + 1, k, n);
            hi = lo;
        } else if (
            bchi == BCType::foextrap || bchi == BCType::hoextrap ||
            bchi == BCType::reflect_even ||
            (bchi == BCType::ext_dir && vedge(i, j, k) > 0)) {
            hi = lo;

        } else if (bchi == BCType::reflect_odd) {
            lo = 0.;
            hi = lo;
        }
    } else {
        return;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_cc_zbc_lo(
    const int i,
    const int j,
    const int k,
    const int n,
    const amrex::Array4<const amrex::Real>& s,
    amrex::Real& lo,
    amrex::Real& hi,
    const amrex::Array4<const amrex::Real>& wedge,
    const int bclo,
    const int domlo)
{
    using namespace amrex;

    if (k == domlo) {
        if (bclo == BCType::ext_dir && wedge(i, j, k) >= 0.) {
            hi = s(i, j, domlo - 1, n);
            lo = hi;
        }

        else if (
            bclo == BCType::foextrap || bclo == BCType::hoextrap ||
            bclo == BCType::reflect_even ||
            (bclo == BCType::ext_dir && wedge(i, j, k) < 0)) {
            lo = hi;

        } else if (bclo == BCType::reflect_odd) {
            hi = 0.;
            lo = hi;
        }
    } else {
        return;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_cc_zbc_hi(
    const int i,
    const int j,
    const int k,
    const int n,
    const amrex::Array4<const amrex::Real>& s,
    amrex::Real& lo,
    amrex::Real& hi,
    const amrex::Array4<const amrex::Real>& wedge,
    const int bchi,
    const int domhi)
{
    using namespace amrex;

    if (k == domhi + 1) {
        if (bchi == BCType::ext_dir && wedge(i, j, k) <= 0.) {
            lo = s(i, j, domhi + 1, n);
            hi = lo;
        } else if (
            bchi == BCType::foextrap || bchi == BCType::hoextrap ||
            bchi == BCType::reflect_even ||
            (bchi == BCType::ext_dir && wedge(i, j, k) > 0)) {
            hi = lo;

        } else if (bchi == BCType::reflect_odd) {
            lo = 0.;
            hi = lo;
        }
    } else {
        return;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_corner_couple_yx(
    amrex::Real& lo1,
    amrex::Real& hi1,
    int i,
    int j,
    int k,
    int n,
    amrex::Real dt,
    amrex::Real dx,
    bool iconserv,
    amrex::Real lo,
    amrex::Real hi,
    amrex::Array4<amrex::Real const> const& s,
    amrex::Array4<amrex::Real const> const& mac,
    amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo -
              dt / (3.e0 * dx) *
                  (state(i + 1, j - 1, k, n) * mac(i + 1, j - 1, k) -
                   state(i, j - 1, k, n) * mac(i, j - 1, k))
              // divu = 0
              //- dt/(3.e0)   *s(i,j-1,k,n)*divu_cc(i,j-1,k)
              + dt / (3.e0 * dx) * s(i, j - 1, k, n) *
                    (mac(i + 1, j - 1, k) - mac(i, j - 1, k));

        hi1 = hi -
              dt / (3.e0 * dx) *
                  (state(i + 1, j, k, n) * mac(i + 1, j, k) -
                   state(i, j, k, n) * mac(i, j, k))
              // divu = 0
              //- dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
              + dt / (3.e0 * dx) * s(i, j, k, n) *
                    (mac(i + 1, j, k) - mac(i, j, k));
    } else {
        lo1 = lo - dt / (6.e0 * dx) *
                       (mac(i + 1, j - 1, k) + mac(i, j - 1, k)) *
                       (state(i + 1, j - 1, k, n) - state(i, j - 1, k, n));
        hi1 = hi - dt / (6.e0 * dx) * (mac(i + 1, j, k) + mac(i, j, k)) *
                       (state(i + 1, j, k, n) - state(i, j, k, n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_corner_couple_zx(
    amrex::Real& lo1,
    amrex::Real& hi1,
    int i,
    int j,
    int k,
    int n,
    amrex::Real dt,
    amrex::Real dx,
    bool iconserv,
    amrex::Real lo,
    amrex::Real hi,
    amrex::Array4<amrex::Real const> const& s,
    amrex::Array4<amrex::Real const> const& mac,
    amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo -
              dt / (3.e0 * dx) *
                  (state(i + 1, j, k - 1, n) * mac(i + 1, j, k - 1) -
                   state(i, j, k - 1, n) * mac(i, j, k - 1))
              // divu = 0
              //- dt/(3.e0)   *s(i,j,k-1,n)*divu_cc(i,j,k-1)
              + dt / (3.e0 * dx) * s(i, j, k - 1, n) *
                    (mac(i + 1, j, k - 1) - mac(i, j, k - 1));

        hi1 = hi -
              dt / (3.e0 * dx) *
                  (state(i + 1, j, k, n) * mac(i + 1, j, k) -
                   state(i, j, k, n) * mac(i, j, k))
              // divu = 0
              //- dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
              + dt / (3.e0 * dx) * s(i, j, k, n) *
                    (mac(i + 1, j, k) - mac(i, j, k));
    } else {
        lo1 = lo - dt / (6.e0 * dx) *
                       (mac(i + 1, j, k - 1) + mac(i, j, k - 1)) *
                       (state(i + 1, j, k - 1, n) - state(i, j, k - 1, n));
        hi1 = hi - dt / (6.e0 * dx) * (mac(i + 1, j, k) + mac(i, j, k)) *
                       (state(i + 1, j, k, n) - state(i, j, k, n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_corner_couple_xy(
    amrex::Real& lo1,
    amrex::Real& hi1,
    int i,
    int j,
    int k,
    int n,
    amrex::Real dt,
    amrex::Real dy,
    bool iconserv,
    amrex::Real lo,
    amrex::Real hi,
    amrex::Array4<amrex::Real const> const& s,
    amrex::Array4<amrex::Real const> const& mac,
    amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo -
              dt / (3.e0 * dy) *
                  (state(i - 1, j + 1, k, n) * mac(i - 1, j + 1, k) -
                   state(i - 1, j, k, n) * mac(i - 1, j, k))
              // divu = 0
              //- dt/3.e0     *s(i-1,j,k,n)*divu_cc(i-1,j,k)
              + dt / (3.e0 * dy) * s(i - 1, j, k, n) *
                    (mac(i - 1, j + 1, k) - mac(i - 1, j, k));

        hi1 = hi -
              dt / (3.e0 * dy) *
                  (state(i, j + 1, k, n) * mac(i, j + 1, k) -
                   state(i, j, k, n) * mac(i, j, k))
              // divu = 0
              //- dt/3.e0     *s(i,j,k,n)*divu_cc(i,j,k)
              + dt / (3.e0 * dy) * s(i, j, k, n) *
                    (mac(i, j + 1, k) - mac(i, j, k));
    } else {
        lo1 = lo - dt / (6.e0 * dy) *
                       (mac(i - 1, j + 1, k) + mac(i - 1, j, k)) *
                       (state(i - 1, j + 1, k, n) - state(i - 1, j, k, n));
        hi1 = hi - dt / (6.e0 * dy) * (mac(i, j + 1, k) + mac(i, j, k)) *
                       (state(i, j + 1, k, n) - state(i, j, k, n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_corner_couple_zy(
    amrex::Real& lo1,
    amrex::Real& hi1,
    int i,
    int j,
    int k,
    int n,
    amrex::Real dt,
    amrex::Real dy,
    bool iconserv,
    amrex::Real lo,
    amrex::Real hi,
    amrex::Array4<amrex::Real const> const& s,
    amrex::Array4<amrex::Real const> const& mac,
    amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo -
              dt / (3.e0 * dy) *
                  (state(i, j + 1, k - 1, n) * mac(i, j + 1, k - 1) -
                   state(i, j, k - 1, n) * mac(i, j, k - 1))
              // divu = 0
              //- dt/(3.e0)   *s(i,j,k-1,n)*divu_cc(i,j,k-1)
              + dt / (3.e0 * dy) * s(i, j, k - 1, n) *
                    (mac(i, j + 1, k - 1) - mac(i, j, k - 1));

        hi1 = hi -
              dt / (3.e0 * dy) *
                  (state(i, j + 1, k, n) * mac(i, j + 1, k) -
                   state(i, j, k, n) * mac(i, j, k))
              // divu = 0
              //- dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
              + dt / (3.e0 * dy) * s(i, j, k, n) *
                    (mac(i, j + 1, k) - mac(i, j, k));
    } else {
        lo1 = lo - dt / (6.e0 * dy) *
                       (mac(i, j + 1, k - 1) + mac(i, j, k - 1)) *
                       (state(i, j + 1, k - 1, n) - state(i, j, k - 1, n));
        hi1 = hi - dt / (6.e0 * dy) * (mac(i, j + 1, k) + mac(i, j, k)) *
                       (state(i, j + 1, k, n) - state(i, j, k, n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_corner_couple_xz(
    amrex::Real& lo1,
    amrex::Real& hi1,
    int i,
    int j,
    int k,
    int n,
    amrex::Real dt,
    amrex::Real dz,
    bool iconserv,
    amrex::Real lo,
    amrex::Real hi,
    amrex::Array4<amrex::Real const> const& s,
    amrex::Array4<amrex::Real const> const& mac,
    amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo -
              dt / (3.e0 * dz) *
                  (state(i - 1, j, k + 1, n) * mac(i - 1, j, k + 1) -
                   state(i - 1, j, k, n) * mac(i - 1, j, k))
              // divu = 0
              //- dt/(3.e0)   *s(i-1,j,k,n)*divu_cc(i-1,j,k)
              + dt / (3.e0 * dz) * s(i - 1, j, k, n) *
                    (mac(i - 1, j, k + 1) - mac(i - 1, j, k));

        hi1 = hi -
              dt / (3.e0 * dz) *
                  (state(i, j, k + 1, n) * mac(i, j, k + 1) -
                   state(i, j, k, n) * mac(i, j, k))
              // divu = 0
              //- dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
              + dt / (3.e0 * dz) * s(i, j, k, n) *
                    (mac(i, j, k + 1) - mac(i, j, k));
    } else {
        lo1 = lo - dt / (6.e0 * dz) *
                       (mac(i - 1, j, k + 1) + mac(i - 1, j, k)) *
                       (state(i - 1, j, k + 1, n) - state(i - 1, j, k, n));
        hi1 = hi - dt / (6.e0 * dz) * (mac(i, j, k + 1) + mac(i, j, k)) *
                       (state(i, j, k + 1, n) - state(i, j, k, n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_corner_couple_yz(
    amrex::Real& lo1,
    amrex::Real& hi1,
    int i,
    int j,
    int k,
    int n,
    amrex::Real dt,
    amrex::Real dz,
    bool iconserv,
    amrex::Real lo,
    amrex::Real hi,
    amrex::Array4<amrex::Real const> const& s,
    amrex::Array4<amrex::Real const> const& mac,
    amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo -
              dt / (3.e0 * dz) *
                  (state(i, j - 1, k + 1, n) * mac(i, j - 1, k + 1) -
                   state(i, j - 1, k, n) * mac(i, j - 1, k))
              // divu = 0
              //- dt/(3.e0)   *s(i,j-1,k,n)*divu_cc(i,j-1,k)
              + dt / (3.e0 * dz) * s(i, j - 1, k, n) *
                    (mac(i, j - 1, k + 1) - mac(i, j - 1, k));

        hi1 = hi -
              dt / (3.e0 * dz) *
                  (state(i, j, k + 1, n) * mac(i, j, k + 1) -
                   state(i, j, k, n) * mac(i, j, k))
              // divu = 0
              //- dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
              + dt / (3.e0 * dz) * s(i, j, k, n) *
                    (mac(i, j, k + 1) - mac(i, j, k));
    } else {
        lo1 = lo - dt / (6.e0 * dz) *
                       (mac(i, j - 1, k + 1) + mac(i, j - 1, k)) *
                       (state(i, j - 1, k + 1, n) - state(i, j - 1, k, n));
        hi1 = hi - dt / (6.e0 * dz) * (mac(i, j, k + 1) + mac(i, j, k)) *
                       (state(i, j, k + 1, n) - state(i, j, k, n));
    }
}

} // namespace
#endif
