#ifndef GODUNOV_PPM_NOLIM_H
#define GODUNOV_PPM_NOLIM_H

#include <AMReX_Gpu.H>
#include <AMReX_BCRec.H>

/* This header file contains the inlined __host__ __device__ functions required
   for the scalar advection routines for 3D Godunov without a limiter. It also
   contains function declarations for controlling host functions.  */

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_ppm_xbc_nolim(
    const int i,
    const int j,
    const int k,
    const int n,
    amrex::Real& sm,
    amrex::Real& sp,
    amrex::Real& sedge1,
    amrex::Real& sedge2,
    const amrex::Array4<const amrex::Real>& s,
    const int bclo,
    const int bchi,
    const int domlo,
    const int domhi)
{
    using namespace amrex;

    if (bclo == BCType::ext_dir || bclo == BCType::hoextrap) {
        if (i == domlo) {
            sedge2 = -0.2 * s(domlo - 1, j, k, n) + 0.75 * s(domlo, j, k, n) +
                     0.5 * s(domlo + 1, j, k, n) - 0.05 * s(domlo + 2, j, k, n);

            sm = s(domlo - 1, j, k, n);
            sp = sedge2;

        } else if (i == domlo + 1) {

            sedge1 = -0.2 * s(domlo - 1, j, k, n) + 0.75 * s(domlo, j, k, n) +
                     0.5 * s(domlo + 1, j, k, n) - 0.05 * s(domlo + 2, j, k, n);

            sp = sedge2;
            sm = sedge1;
        }
    }

    if (bchi == BCType::ext_dir || bchi == BCType::hoextrap) {
        if (i == domhi) {
            sedge1 = -0.2 * s(domhi + 1, j, k, n) + 0.75 * s(domhi, j, k, n) +
                     0.5 * s(domhi - 1, j, k, n) - 0.05 * s(domhi - 2, j, k, n);

            sp = s(domhi + 1, j, k, n);
            sm = sedge1;

        } else if (i == domhi - 1) {

            sedge2 = -0.2 * s(domhi + 1, j, k, n) + 0.75 * s(domhi, j, k, n) +
                     0.5 * s(domhi - 1, j, k, n) - 0.05 * s(domhi - 2, j, k, n);

            sp = sedge2;
            sm = sedge1;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_ppm_ybc_nolim(
    const int i,
    const int j,
    const int k,
    const int n,
    amrex::Real& sm,
    amrex::Real& sp,
    amrex::Real& sedge1,
    amrex::Real& sedge2,
    const amrex::Array4<const amrex::Real>& s,
    const int bclo,
    const int bchi,
    const int domlo,
    const int domhi)
{
    using namespace amrex;

    if (bclo == BCType::ext_dir || bclo == BCType::hoextrap) {
        if (j == domlo) {
            sedge2 = -0.2 * s(i, domlo - 1, k, n) + 0.75 * s(i, domlo, k, n) +
                     0.5 * s(i, domlo + 1, k, n) - 0.05 * s(i, domlo + 2, k, n);

            sm = s(i, domlo - 1, k, n);
            sp = sedge2;

        } else if (j == domlo + 1) {

            sedge1 = -0.2 * s(i, domlo - 1, k, n) + 0.75 * s(i, domlo, k, n) +
                     0.5 * s(i, domlo + 1, k, n) - 0.05 * s(i, domlo + 2, k, n);

            sp = sedge2;
            sm = sedge1;
        }
    }

    if (bchi == BCType::ext_dir || bchi == BCType::hoextrap) {
        if (j == domhi) {
            sedge1 = -0.2 * s(i, domhi + 1, k, n) + 0.75 * s(i, domhi, k, n) +
                     0.5 * s(i, domhi - 1, k, n) - 0.05 * s(i, domhi - 2, k, n);

            sp = s(i, domhi + 1, k, n);
            sm = sedge1;

        } else if (j == domhi - 1) {

            sedge2 = -0.2 * s(i, domhi + 1, k, n) + 0.75 * s(i, domhi, k, n) +
                     0.5 * s(i, domhi - 1, k, n) - 0.05 * s(i, domhi - 2, k, n);

            sp = sedge2;
            sm = sedge1;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_ppm_zbc_nolim(
    const int i,
    const int j,
    const int k,
    const int n,
    amrex::Real& sm,
    amrex::Real& sp,
    amrex::Real& sedge1,
    amrex::Real& sedge2,
    const amrex::Array4<const amrex::Real>& s,
    const int bclo,
    const int bchi,
    const int domlo,
    const int domhi)
{
    using namespace amrex;

    if (bclo == BCType::ext_dir || bclo == BCType::hoextrap) {

        if (k == domlo) {
            sedge2 = -0.2 * s(i, j, domlo - 1, n) + 0.75 * s(i, j, domlo, n) +
                     0.5 * s(i, j, domlo + 1, n) - 0.05 * s(i, j, domlo + 2, n);

            sm = s(i, j, domlo - 1, n);
            sp = sedge2;

        } else if (k == domlo + 1) {

            sedge1 = -0.2 * s(i, j, domlo - 1, n) + 0.75 * s(i, j, domlo, n) +
                     0.5 * s(i, j, domlo + 1, n) - 0.05 * s(i, j, domlo + 2, n);

            sp = sedge2;
            sm = sedge1;
        }
    }

    if (bchi == BCType::ext_dir || bchi == BCType::hoextrap) {
        if (k == domhi) {
            sedge1 = -0.2 * s(i, j, domhi + 1, n) + 0.75 * s(i, j, domhi, n) +
                     0.5 * s(i, j, domhi - 1, n) - 0.05 * s(i, j, domhi - 2, n);

            sp = s(i, j, domhi + 1, n);
            sm = sedge1;

        } else if (k == domhi - 1) {

            sedge2 = -0.2 * s(i, j, domhi + 1, n) + 0.75 * s(i, j, domhi, n) +
                     0.5 * s(i, j, domhi - 1, n) - 0.05 * s(i, j, domhi - 2, n);

            sp = sedge2;
            sm = sedge1;
        }
    }
}

// Right now only ppm type 1 is supported on GPU
// This version is called before the MAC projection, when we use the
// cell-centered velocity
//      for upwinding
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_ppm_pred_x_nolim(
    const int i,
    const int j,
    const int k,
    const int n,
    const amrex::Real dtdx,
    const amrex::Real v_ad,
    const amrex::Array4<const amrex::Real>& S,
    const amrex::Array4<amrex::Real>& Im,
    const amrex::Array4<amrex::Real>& Ip,
    const amrex::BCRec& bc,
    const int domlo,
    const int domhi)
{
    using namespace amrex;

    constexpr amrex::Real small_vel = 1e-10;
    constexpr amrex::Real sixth = 1.0 / 6.0;

    amrex::Real sedge1, sedge2, sp, sm;

    amrex::Real sm2 = S(i - 2, j, k, n);
    amrex::Real sm1 = S(i - 1, j, k, n);
    amrex::Real s0 = S(i, j, k, n);
    amrex::Real sp1 = S(i + 1, j, k, n);
    amrex::Real sp2 = S(i + 2, j, k, n);

    amrex::Real d1 = 0.5 * (sp1 - sm1);
    amrex::Real d2 = 0.5 * (s0 - sm2);

    sedge1 = 0.5e0 * (s0 + sm1) - sixth * (d1 - d2);

    d1 = 0.5 * (sp2 - s0);
    d2 = 0.5 * (sp1 - sm1);

    sedge2 = 0.5e0 * (sp1 + s0) - sixth * (d1 - d2);

    sm = sedge1;
    sp = sedge2;

    Godunov_ppm_xbc_nolim(
        i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(0), bc.hi(0), domlo,
        domhi);

    amrex::Real s6 = 6.0 * s0 - 3.0 * (sm + sp);

    amrex::Real sigma = std::abs(v_ad) * dtdx;

    // upwind
    if (v_ad > small_vel) {
        Ip(i, j, k, n) =
            sp - (0.5 * sigma) * ((sp - sm) - (1.0 - 2.0 / 3.0 * sigma) * s6);
        Im(i, j, k, n) = S(i, j, k, n);
    } else if (v_ad < -small_vel) {
        Ip(i, j, k, n) = S(i, j, k, n);
        Im(i, j, k, n) =
            sm + (0.5 * sigma) * ((sp - sm) + (1.0 - 2.0 / 3.0 * sigma) * s6);
    } else {
        Ip(i, j, k, n) = S(i, j, k, n);
        Im(i, j, k, n) = S(i, j, k, n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_ppm_pred_y_nolim(
    const int i,
    const int j,
    const int k,
    const int n,
    const amrex::Real dtdy,
    const amrex::Real v_ad,
    const amrex::Array4<const amrex::Real>& S,
    const amrex::Array4<amrex::Real>& Im,
    const amrex::Array4<amrex::Real>& Ip,
    const amrex::BCRec& bc,
    const int domlo,
    const int domhi)
{
    using namespace amrex;

    constexpr amrex::Real small_vel = 1e-10;
    constexpr amrex::Real sixth = 1.0 / 6.0;

    amrex::Real sedge1, sedge2, sp, sm;

    amrex::Real sm2 = S(i, j - 2, k, n);
    amrex::Real sm1 = S(i, j - 1, k, n);
    amrex::Real s0 = S(i, j, k, n);
    amrex::Real sp1 = S(i, j + 1, k, n);
    amrex::Real sp2 = S(i, j + 2, k, n);

    amrex::Real d1 = 0.5 * (sp1 - sm1);
    amrex::Real d2 = 0.5 * (s0 - sm2);

    sedge1 = 0.5e0 * (s0 + sm1) - sixth * (d1 - d2);

    d1 = 0.5 * (sp2 - s0);
    d2 = 0.5 * (sp1 - sm1);

    sedge2 = 0.5e0 * (sp1 + s0) - sixth * (d1 - d2);

    sm = sedge1;
    sp = sedge2;

    Godunov_ppm_ybc_nolim(
        i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(1), bc.hi(1), domlo,
        domhi);

    amrex::Real s6 = 6.0 * s0 - 3.0 * (sm + sp);

    amrex::Real sigma = std::abs(v_ad) * dtdy;

    if (v_ad > small_vel) {
        Ip(i, j, k, n) =
            sp - (0.5 * sigma) * ((sp - sm) - (1.0 - 2.0 / 3.0 * sigma) * s6);
        Im(i, j, k, n) = S(i, j, k, n);
    } else if (v_ad < -small_vel) {
        Ip(i, j, k, n) = S(i, j, k, n);
        Im(i, j, k, n) =
            sm + (0.5 * sigma) * ((sp - sm) + (1.0 - 2.0 / 3.0 * sigma) * s6);
    } else {
        Ip(i, j, k, n) = S(i, j, k, n);
        Im(i, j, k, n) = S(i, j, k, n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_ppm_pred_z_nolim(
    const int i,
    const int j,
    const int k,
    const int n,
    const amrex::Real dtdz,
    const amrex::Real v_ad,
    const amrex::Array4<const amrex::Real>& S,
    const amrex::Array4<amrex::Real>& Im,
    const amrex::Array4<amrex::Real>& Ip,
    const amrex::BCRec& bc,
    const int domlo,
    const int domhi)
{
    using namespace amrex;

    constexpr amrex::Real small_vel = 1e-10;
    constexpr amrex::Real sixth = 1.0 / 6.0;

    amrex::Real sedge1, sedge2, sp, sm;

    amrex::Real sm2 = S(i, j, k - 2, n);
    amrex::Real sm1 = S(i, j, k - 1, n);
    amrex::Real s0 = S(i, j, k, n);
    amrex::Real sp1 = S(i, j, k + 1, n);
    amrex::Real sp2 = S(i, j, k + 2, n);

    amrex::Real d1 = 0.5 * (sp1 - sm1);
    amrex::Real d2 = 0.5 * (s0 - sm2);

    sedge1 = 0.5e0 * (s0 + sm1) - sixth * (d1 - d2);

    d1 = 0.5 * (sp2 - s0);
    d2 = 0.5 * (sp1 - sm1);

    sedge2 = 0.5e0 * (sp1 + s0) - sixth * (d1 - d2);

    sm = sedge1;
    sp = sedge2;

    Godunov_ppm_zbc_nolim(
        i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(2), bc.hi(2), domlo,
        domhi);

    amrex::Real s6 = 6.0 * s0 - 3.0 * (sm + sp);

    amrex::Real sigma = std::abs(v_ad) * dtdz;

    if (v_ad > small_vel) {
        Ip(i, j, k, n) =
            sp - (0.5 * sigma) * ((sp - sm) - (1.0 - 2.0 / 3.0 * sigma) * s6);
        Im(i, j, k, n) = S(i, j, k, n);
    } else if (v_ad < -small_vel) {
        Ip(i, j, k, n) = S(i, j, k, n);
        Im(i, j, k, n) =
            sm + (0.5 * sigma) * ((sp - sm) + (1.0 - 2.0 / 3.0 * sigma) * s6);
    } else {
        Ip(i, j, k, n) = S(i, j, k, n);
        Im(i, j, k, n) = S(i, j, k, n);
    }
}

// Right now only ppm type 1 is supported on GPU
// This version is called after the MAC projection, when we use the
// MAC-projected velocity
//      for upwinding
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_ppm_fpu_x_nolim(
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
    {
        using namespace amrex;

        constexpr Real small_vel = 1e-10;
        constexpr Real sixth = 1.0 / 6.0;

        amrex::Real sedge1, sedge2, sp, sm;

        amrex::Real sm2 = S(i - 2, j, k, n);
        amrex::Real sm1 = S(i - 1, j, k, n);
        amrex::Real s0 = S(i, j, k, n);
        amrex::Real sp1 = S(i + 1, j, k, n);
        amrex::Real sp2 = S(i + 2, j, k, n);

        amrex::Real d1 = 0.5 * (sp1 - sm1);
        amrex::Real d2 = 0.5 * (s0 - sm2);

        sedge1 = 0.5e0 * (s0 + sm1) - sixth * (d1 - d2);

        d1 = 0.5 * (sp2 - s0);
        d2 = 0.5 * (sp1 - sm1);

        sedge2 = 0.5e0 * (sp1 + s0) - sixth * (d1 - d2);

        sm = sedge1;
        sp = sedge2;

        Godunov_ppm_xbc_nolim(
            i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(0), bc.hi(0), domlo,
            domhi);

        Real s6 = 6.0 * s0 - 3.0 * (sm + sp);

        Real sigmap = std::abs(vel_edge(i + 1, j, k)) * dt / dx;
        Real sigmam = std::abs(vel_edge(i, j, k)) * dt / dx;

        if (vel_edge(i + 1, j, k) > small_vel) {
            Ip = sp - (0.5 * sigmap) *
                          ((sp - sm) - (1.e0 - 2.e0 / 3.e0 * sigmap) * s6);
        } else {
            Ip = S(i, j, k, n);
        }

        if (vel_edge(i, j, k) < -small_vel) {
            Im = sm + (0.5 * sigmam) *
                          ((sp - sm) + (1.e0 - 2.e0 / 3.e0 * sigmam) * s6);
        } else {
            Im = S(i, j, k, n);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_ppm_fpu_y_nolim(
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
    {
        using namespace amrex;

        constexpr amrex::Real small_vel = 1e-10;
        constexpr amrex::Real sixth = 1.0 / 6.0;

        Real sedge1, sedge2, sp, sm;

        Real sm2 = S(i, j - 2, k, n);
        Real sm1 = S(i, j - 1, k, n);
        Real s0 = S(i, j, k, n);
        Real sp1 = S(i, j + 1, k, n);
        Real sp2 = S(i, j + 2, k, n);

        Real d1 = 0.5 * (sp1 - sm1);
        Real d2 = 0.5 * (s0 - sm2);

        sedge1 = 0.5e0 * (s0 + sm1) - sixth * (d1 - d2);

        d1 = 0.5 * (sp2 - s0);
        d2 = 0.5 * (sp1 - sm1);

        sedge2 = 0.5e0 * (sp1 + s0) - sixth * (d1 - d2);

        sm = sedge1;
        sp = sedge2;

        Godunov_ppm_ybc_nolim(
            i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(1), bc.hi(1), domlo,
            domhi);

        amrex::Real s6 = 6.0 * s0 - 3.0 * (sm + sp);

        Real sigmap = std::abs(vel_edge(i, j + 1, k)) * dt / dx;
        Real sigmam = std::abs(vel_edge(i, j, k)) * dt / dx;

        if (vel_edge(i, j + 1, k) > small_vel) {
            Ip = sp - (0.5 * sigmap) *
                          ((sp - sm) - (1.e0 - 2.e0 / 3.e0 * sigmap) * s6);
        } else {
            Ip = S(i, j, k, n);
        }

        if (vel_edge(i, j, k) < -small_vel) {
            Im = sm + (0.5 * sigmam) *
                          ((sp - sm) + (1.e0 - 2.e0 / 3.e0 * sigmam) * s6);
        } else {
            Im = S(i, j, k, n);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_ppm_fpu_z_nolim(
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
    {
        using namespace amrex;

        constexpr amrex::Real small_vel = 1e-10;
        constexpr amrex::Real sixth = 1.0 / 6.0;

        Real sedge1, sedge2, sp, sm;

        Real sm2 = S(i, j, k - 2, n);
        Real sm1 = S(i, j, k - 1, n);
        Real s0 = S(i, j, k, n);
        Real sp1 = S(i, j, k + 1, n);
        Real sp2 = S(i, j, k + 2, n);

        Real d1 = 0.5 * (sp1 - sm1);
        Real d2 = 0.5 * (s0 - sm2);

        sedge1 = 0.5e0 * (s0 + sm1) - sixth * (d1 - d2);

        d1 = 0.5 * (sp2 - s0);
        d2 = 0.5 * (sp1 - sm1);

        sedge2 = 0.5e0 * (sp1 + s0) - sixth * (d1 - d2);

        sm = sedge1;
        sp = sedge2;

        Godunov_ppm_zbc_nolim(
            i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(2), bc.hi(2), domlo,
            domhi);

        Real s6 = 6.0 * s0 - 3.0 * (sm + sp);
        Real sigmap = std::abs(vel_edge(i, j, k + 1)) * dt / dx;
        Real sigmam = std::abs(vel_edge(i, j, k)) * dt / dx;

        if (vel_edge(i, j, k + 1) > small_vel) {
            Ip = sp - (0.5 * sigmap) *
                          ((sp - sm) - (1.e0 - 2.e0 / 3.e0 * sigmap) * s6);
        } else {
            Ip = S(i, j, k, n);
        }

        if (vel_edge(i, j, k) < -small_vel) {
            Im = sm + (0.5 * sigmam) *
                          ((sp - sm) + (1.e0 - 2.e0 / 3.e0 * sigmam) * s6);
        } else {
            Im = S(i, j, k, n);
        }
    }
}

} // namespace
#endif
