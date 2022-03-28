#ifndef GODUNOV_MINMOD_H
#define GODUNOV_MINMOD_H

#include <AMReX_Gpu.H>
#include <AMReX_BCRec.H>

/* This header file contains the inlined __host__ __device__ functions required
   for the scalar advection routines for 3D Godunov. It also contains function
   declarations for controlling host functions.  */

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_minmod_xbc(
    const int i,
    const int j,
    const int k,
    const int n,
    amrex::Real& sm,
    amrex::Real& sp,
    const amrex::Real& sedge1,
    const amrex::Real& sedge2,
    const amrex::Array4<const amrex::Real>& s,
    const int bclo,
    const int bchi,
    const int domlo,
    const int domhi)
{
    using namespace amrex;

    if (bclo == BCType::ext_dir || bclo == BCType::hoextrap) {
        if (i == domlo) {
            sp = -0.2 * s(domlo - 1, j, k, n) + 0.75 * s(domlo, j, k, n) +
                 0.5 * s(domlo + 1, j, k, n) - 0.05 * s(domlo + 2, j, k, n);

            sm = s(domlo - 1, j, k, n);

        } else if (i == domlo + 1) {

            sm = -0.2 * s(domlo - 1, j, k, n) + 0.75 * s(domlo, j, k, n) +
                 0.5 * s(domlo + 1, j, k, n) - 0.05 * s(domlo + 2, j, k, n);

            sp = sedge2;
        }
    }

    if (bchi == BCType::ext_dir || bchi == BCType::hoextrap) {
        if (i == domhi) {
            sm = -0.2 * s(domhi + 1, j, k, n) + 0.75 * s(domhi, j, k, n) +
                 0.5 * s(domhi - 1, j, k, n) - 0.05 * s(domhi - 2, j, k, n);

            sp = s(domhi + 1, j, k, n);

        } else if (i == domhi - 1) {

            sm = -0.2 * s(domhi + 1, j, k, n) + 0.75 * s(domhi, j, k, n) +
                 0.5 * s(domhi - 1, j, k, n) - 0.05 * s(domhi - 2, j, k, n);

            sp = sedge1;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_minmod_ybc(
    const int i,
    const int j,
    const int k,
    const int n,
    amrex::Real& sm,
    amrex::Real& sp,
    const amrex::Real& sedge1,
    const amrex::Real& sedge2,
    const amrex::Array4<const amrex::Real>& s,
    const int bclo,
    const int bchi,
    const int domlo,
    const int domhi)
{
    using namespace amrex;

    if (bclo == BCType::ext_dir || bclo == BCType::hoextrap) {
        if (j == domlo) {
            sp = -0.2 * s(i, domlo - 1, k, n) + 0.75 * s(i, domlo, k, n) +
                 0.5 * s(i, domlo + 1, k, n) - 0.05 * s(i, domlo + 2, k, n);

            sm = s(i, domlo - 1, k, n);

        } else if (j == domlo + 1) {

            sm = -0.2 * s(i, domlo - 1, k, n) + 0.75 * s(i, domlo, k, n) +
                 0.5 * s(i, domlo + 1, k, n) - 0.05 * s(i, domlo + 2, k, n);

            sp = sedge2;
        }
    }

    if (bchi == BCType::ext_dir || bchi == BCType::hoextrap) {
        if (j == domhi) {
            sm = -0.2 * s(i, domhi + 1, k, n) + 0.75 * s(i, domhi, k, n) +
                 0.5 * s(i, domhi - 1, k, n) - 0.05 * s(i, domhi - 2, k, n);

            sp = s(i, domhi + 1, k, n);

        } else if (j == domhi - 1) {

            sp = -0.2 * s(i, domhi + 1, k, n) + 0.75 * s(i, domhi, k, n) +
                 0.5 * s(i, domhi - 1, k, n) - 0.05 * s(i, domhi - 2, k, n);

            sm = sedge1;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_minmod_zbc(
    const int i,
    const int j,
    const int k,
    const int n,
    amrex::Real& sm,
    amrex::Real& sp,
    const amrex::Real& sedge1,
    const amrex::Real& sedge2,
    const amrex::Array4<const amrex::Real>& s,
    const int bclo,
    const int bchi,
    const int domlo,
    const int domhi)
{
    using namespace amrex;

    if (bclo == BCType::ext_dir || bclo == BCType::hoextrap) {

        if (k == domlo) {
            sp = -0.2 * s(i, j, domlo - 1, n) + 0.75 * s(i, j, domlo, n) +
                 0.5 * s(i, j, domlo + 1, n) - 0.05 * s(i, j, domlo + 2, n);

            sm = s(i, j, domlo - 1, n);

        } else if (k == domlo + 1) {

            sm = -0.2 * s(i, j, domlo - 1, n) + 0.75 * s(i, j, domlo, n) +
                 0.5 * s(i, j, domlo + 1, n) - 0.05 * s(i, j, domlo + 2, n);

            sp = sedge2;
        }
    }

    if (bchi == BCType::ext_dir || bchi == BCType::hoextrap) {
        if (k == domhi) {
            sm = -0.2 * s(i, j, domhi + 1, n) + 0.75 * s(i, j, domhi, n) +
                 0.5 * s(i, j, domhi - 1, n) - 0.05 * s(i, j, domhi - 2, n);

            sp = s(i, j, domhi + 1, n);

        } else if (k == domhi - 1) {

            sp = -0.2 * s(i, j, domhi + 1, n) + 0.75 * s(i, j, domhi, n) +
                 0.5 * s(i, j, domhi - 1, n) - 0.05 * s(i, j, domhi - 2, n);

            sm = sedge1;
        }
    }
}

// This version is called after the MAC projection, when we use the
// MAC-projected velocity
//      for upwinding
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_minmod_fpu_x(
    const int i,
    const int j,
    const int k,
    const int n,
    const amrex::Real dt,
    const amrex::Real dx,
    amrex::Real& Im,
    amrex::Real& Ip,
    const amrex::Array4<const amrex::Real>& S,
    const amrex::Array4<const amrex::Real>& vel_edge,
    const amrex::BCRec& bc,
    const int domlo,
    const int domhi)
{

    using namespace amrex;

    constexpr amrex::Real small_vel = 1e-10;

    amrex::Real sm1 = S(i - 1, j, k, n);
    amrex::Real s0 = S(i, j, k, n);
    amrex::Real sp1 = S(i + 1, j, k, n);

    amrex::Real dsp = sp1 - s0;
    amrex::Real dsm = s0 - sm1;

    amrex::Real ds = 0.0;
    if (!(dsp*dsm < 0.0)) {
        if (amrex::Math::abs(dsp) < amrex::Math::abs(dsm)) {
            ds = dsp;
        } else {
            ds = dsm;
        }
    }

    amrex::Real sigmap = amrex::Math::abs(vel_edge(i + 1, j, k)) * dt / dx;
    amrex::Real sigmam = amrex::Math::abs(vel_edge(i, j, k)) * dt / dx;

    if (vel_edge(i + 1, j, k) > small_vel) {
        Ip = s0 + 0.5 * (1.0 - sigmap) * ds;
    } else {
        Ip = S(i, j, k, n);
    }

    if (vel_edge(i, j, k) < -small_vel) {
        Im = s0 - 0.5 * (1.0 - sigmam) * ds;
    } else {
        Im = S(i, j, k, n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_minmod_fpu_y(
    const int i,
    const int j,
    const int k,
    const int n,
    const amrex::Real dt,
    const amrex::Real dx,
    amrex::Real& Im,
    amrex::Real& Ip,
    const amrex::Array4<const amrex::Real>& S,
    const amrex::Array4<const amrex::Real>& vel_edge,
    const amrex::BCRec& bc,
    const int domlo,
    const int domhi)
{

    using namespace amrex;

    constexpr amrex::Real small_vel = 1e-10;

    amrex::Real sm1 = S(i, j - 1, k, n);
    amrex::Real s0 = S(i, j, k, n);
    amrex::Real sp1 = S(i, j + 1, k, n);

    amrex::Real dsp = sp1 - s0;
    amrex::Real dsm = s0 - sm1;

    amrex::Real ds = 0.0;
    if (!(dsp*dsm < 0.0)) {
        if (amrex::Math::abs(dsp) < amrex::Math::abs(dsm)) {
            ds = dsp;
        } else {
            ds = dsm;
        }
    }

    amrex::Real sigmap = amrex::Math::abs(vel_edge(i, j + 1, k)) * dt / dx;
    amrex::Real sigmam = amrex::Math::abs(vel_edge(i, j, k)) * dt / dx;

    if (vel_edge(i, j + 1, k) > small_vel) {
        Ip = s0 + 0.5 * (1.0 - sigmap) * ds;
    } else {
        Ip = S(i, j, k, n);
    }

    if (vel_edge(i, j, k) < -small_vel) {
        Im = s0 - 0.5 * (1.0 - sigmam) * ds;
    } else {
        Im = S(i, j, k, n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_minmod_fpu_z(
    const int i,
    const int j,
    const int k,
    const int n,
    const amrex::Real dt,
    const amrex::Real dx,
    amrex::Real& Im,
    amrex::Real& Ip,
    const amrex::Array4<const amrex::Real>& S,
    const amrex::Array4<const amrex::Real>& vel_edge,
    const amrex::BCRec& bc,
    const int domlo,
    const int domhi)
{

    using namespace amrex;

    constexpr amrex::Real small_vel = 1e-10;

    amrex::Real sm1 = S(i, j, k - 1, n);
    amrex::Real s0 = S(i, j, k, n);
    amrex::Real sp1 = S(i, j, k + 1, n);

    amrex::Real dsp = sp1 - s0;
    amrex::Real dsm = s0 - sm1;

    amrex::Real ds = 0.0;
    if (!(dsp*dsm < 0.0)) {
        if (amrex::Math::abs(dsp) < amrex::Math::abs(dsm)) {
            ds = dsp;
        } else {
            ds = dsm;
        }
    }

    amrex::Real sigmap = amrex::Math::abs(vel_edge(i, j, k + 1)) * dt / dx;
    amrex::Real sigmam = amrex::Math::abs(vel_edge(i, j, k)) * dt / dx;

    if (vel_edge(i, j, k + 1) > small_vel) {
        Ip = s0 + 0.5 * (1.0 - sigmap) * ds;
    } else {
        Ip = S(i, j, k, n);
    }

    if (vel_edge(i, j, k) < -small_vel) {
        Im = s0 - 0.5 * (1.0 - sigmam) * ds;
    } else {
        Im = S(i, j, k, n);
    }
}

} // namespace

#endif
