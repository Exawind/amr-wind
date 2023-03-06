#ifndef GODUNOV_WENO_H
#define GODUNOV_WENO_H

#include <AMReX_Gpu.H>
#include <AMReX_BCRec.H>

/* This header file contains the inlined __host__ __device__ functions required
   for the scalar advection routines for 3D Godunov. It also contains function
   declarations for controlling host functions.  */

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real weno5(
    const amrex::Real sm2,
    const amrex::Real sm1,
    const amrex::Real s,
    const amrex::Real sp1,
    const amrex::Real sp2,
    bool weno_js)
{
    constexpr amrex::Real eps = 1.e-6;

    const amrex::Real beta1 =
        13.0 / 12.0 * (sm2 - 2.0 * sm1 + s) * (sm2 - 2.0 * sm1 + s) +
        0.25 * (sm2 - 4.0 * sm1 + 3.0 * s) * (sm2 - 4.0 * sm1 + 3.0 * s);
    const amrex::Real beta2 =
        13.0 / 12.0 * (sm1 - 2.0 * s + sp1) * (sm1 - 2.0 * s + sp1) +
        0.25 * (sm1 - sp1) * (sm1 - sp1);
    const amrex::Real beta3 =
        13.0 / 12.0 * (s - 2.0 * sp1 + sp2) * (s - 2.0 * sp1 + sp2) +
        0.25 * (3.0 * s - 4.0 * sp1 + sp2) * (3.0 * s - 4.0 * sp1 + sp2);

    amrex::Real omega1, omega2, omega3;

    if (weno_js) {
        omega1 = 0.1 / (eps + beta1);
        omega2 = 0.6 / (eps + beta2);
        omega3 = 0.3 / (eps + beta3);
    } else {
        const amrex::Real t5 = std::abs(beta3 - beta1);
        omega1 = 0.1 * (1.0 + t5 / (eps + beta1));
        omega2 = 0.6 * (1.0 + t5 / (eps + beta2));
        omega3 = 0.3 * (1.0 + t5 / (eps + beta3));
    }

    const amrex::Real omega = omega1 + omega2 + omega3;

    const amrex::Real v_1 = 2.0 * sm2 - 7.0 * sm1 + 11.0 * s;
    const amrex::Real v_2 = -sm1 + 5.0 * s + 2.0 * sp1;
    const amrex::Real v_3 = 2.0 * s + 5.0 * sp1 - sp2;

    return (omega1 * v_1 + omega2 * v_2 + omega3 * v_3) / (6.0 * omega);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_weno_xbc(
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_weno_ybc(
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_weno_zbc(
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

// This version is called before the MAC projection, when we use the
// cell-centered velocity
//      for upwinding
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_weno_pred_x(
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
    const int domhi,
    const bool weno_js)
{

    using namespace amrex;

    constexpr amrex::Real small_vel = 1e-10;

    amrex::Real sm2 = S(i - 2, j, k, n);
    amrex::Real sm1 = S(i - 1, j, k, n);
    amrex::Real s0 = S(i, j, k, n);
    amrex::Real sp1 = S(i + 1, j, k, n);
    amrex::Real sp2 = S(i + 2, j, k, n);

    // right of i-1/2
    amrex::Real sedge1 = weno5(sp2, sp1, s0, sm1, sm2, weno_js); // NOLINT
    // left of i+1/2
    amrex::Real sedge2 = weno5(sm2, sm1, s0, sp1, sp2, weno_js); // NOLINT

    amrex::Real sm = sedge1;
    amrex::Real sp = sedge2;

    Godunov_weno_xbc(
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_weno_pred_y(
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
    const int domhi,
    const bool weno_js)
{
    using namespace amrex;

    constexpr amrex::Real small_vel = 1e-10;

    amrex::Real sm2 = S(i, j - 2, k, n);
    amrex::Real sm1 = S(i, j - 1, k, n);
    amrex::Real s0 = S(i, j, k, n);
    amrex::Real sp1 = S(i, j + 1, k, n);
    amrex::Real sp2 = S(i, j + 2, k, n);

    // right of j-1/2
    amrex::Real sedge1 = weno5(sp2, sp1, s0, sm1, sm2, weno_js); // NOLINT
    // left of j+1/2
    amrex::Real sedge2 = weno5(sm2, sm1, s0, sp1, sp2, weno_js); // NOLINT

    amrex::Real sm = sedge1;
    amrex::Real sp = sedge2;

    Godunov_weno_ybc(
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_weno_pred_z(
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
    const int domhi,
    const bool weno_js)
{
    using namespace amrex;

    constexpr amrex::Real small_vel = 1e-10;

    amrex::Real sm2 = S(i, j, k - 2, n);
    amrex::Real sm1 = S(i, j, k - 1, n);
    amrex::Real s0 = S(i, j, k, n);
    amrex::Real sp1 = S(i, j, k + 1, n);
    amrex::Real sp2 = S(i, j, k + 2, n);

    // right of k-1/2
    amrex::Real sedge1 = weno5(sp2, sp1, s0, sm1, sm2, weno_js); // NOLINT
    // left of k+1/2
    amrex::Real sedge2 = weno5(sm2, sm1, s0, sp1, sp2, weno_js); // NOLINT

    amrex::Real sm = sedge1;
    amrex::Real sp = sedge2;

    Godunov_weno_zbc(
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

// This version is called after the MAC projection, when we use the
// MAC-projected velocity
//      for upwinding
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_weno_fpu_x(
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
    const int domhi,
    const bool weno_js)
{

    using namespace amrex;

    constexpr amrex::Real small_vel = 1e-10;

    amrex::Real sm2 = S(i - 2, j, k, n);
    amrex::Real sm1 = S(i - 1, j, k, n);
    amrex::Real s0 = S(i, j, k, n);
    amrex::Real sp1 = S(i + 1, j, k, n);
    amrex::Real sp2 = S(i + 2, j, k, n);

    // right of i-1/2
    amrex::Real sedge1 = weno5(sp2, sp1, s0, sm1, sm2, weno_js); // NOLINT
    // left of i+1/2
    amrex::Real sedge2 = weno5(sm2, sm1, s0, sp1, sp2, weno_js); // NOLINT

    amrex::Real sm = sedge1;
    amrex::Real sp = sedge2;

    Godunov_weno_xbc(
        i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(0), bc.hi(0), domlo,
        domhi);

    amrex::Real s6 = 6.0 * s0 - 3.0 * (sm + sp);

    amrex::Real sigmap = std::abs(vel_edge(i + 1, j, k)) * dt / dx;
    amrex::Real sigmam = std::abs(vel_edge(i, j, k)) * dt / dx;

    if (vel_edge(i + 1, j, k) > small_vel) {
        Ip = sp -
             (0.5 * sigmap) * ((sp - sm) - (1.e0 - 2.e0 / 3.e0 * sigmap) * s6);
    } else {
        Ip = S(i, j, k, n);
    }

    if (vel_edge(i, j, k) < -small_vel) {
        Im = sm +
             (0.5 * sigmam) * ((sp - sm) + (1.e0 - 2.e0 / 3.e0 * sigmam) * s6);
    } else {
        Im = S(i, j, k, n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_weno_fpu_y(
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
    const int domhi,
    const bool weno_js)
{

    using namespace amrex;

    constexpr amrex::Real small_vel = 1e-10;

    amrex::Real sm2 = S(i, j - 2, k, n);
    amrex::Real sm1 = S(i, j - 1, k, n);
    amrex::Real s0 = S(i, j, k, n);
    amrex::Real sp1 = S(i, j + 1, k, n);
    amrex::Real sp2 = S(i, j + 2, k, n);

    // right of j-1/2
    amrex::Real sedge1 = weno5(sp2, sp1, s0, sm1, sm2, weno_js); // NOLINT
    // left of j+1/2
    amrex::Real sedge2 = weno5(sm2, sm1, s0, sp1, sp2, weno_js); // NOLINT

    amrex::Real sm = sedge1;
    amrex::Real sp = sedge2;

    Godunov_weno_ybc(
        i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(1), bc.hi(1), domlo,
        domhi);

    amrex::Real s6 = 6.0 * s0 - 3.0 * (sm + sp);

    amrex::Real sigmap = std::abs(vel_edge(i, j + 1, k)) * dt / dx;
    amrex::Real sigmam = std::abs(vel_edge(i, j, k)) * dt / dx;

    if (vel_edge(i, j + 1, k) > small_vel) {
        Ip = sp -
             (0.5 * sigmap) * ((sp - sm) - (1.e0 - 2.e0 / 3.e0 * sigmap) * s6);
    } else {
        Ip = S(i, j, k, n);
    }

    if (vel_edge(i, j, k) < -small_vel) {
        Im = sm +
             (0.5 * sigmam) * ((sp - sm) + (1.e0 - 2.e0 / 3.e0 * sigmam) * s6);
    } else {
        Im = S(i, j, k, n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_weno_fpu_z(
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
    const int domhi,
    const bool weno_js)
{

    using namespace amrex;

    constexpr amrex::Real small_vel = 1e-10;

    amrex::Real sm2 = S(i, j, k - 2, n);
    amrex::Real sm1 = S(i, j, k - 1, n);
    amrex::Real s0 = S(i, j, k, n);
    amrex::Real sp1 = S(i, j, k + 1, n);
    amrex::Real sp2 = S(i, j, k + 2, n);

    // right of k-1/2
    amrex::Real sedge1 = weno5(sp2, sp1, s0, sm1, sm2, weno_js); // NOLINT
    // left of k+1/2
    amrex::Real sedge2 = weno5(sm2, sm1, s0, sp1, sp2, weno_js); // NOLINT

    amrex::Real sm = sedge1;
    amrex::Real sp = sedge2;

    Godunov_weno_zbc(
        i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(2), bc.hi(2), domlo,
        domhi);

    amrex::Real s6 = 6.0 * s0 - 3.0 * (sm + sp);

    amrex::Real sigmap = std::abs(vel_edge(i, j, k + 1)) * dt / dx;
    amrex::Real sigmam = std::abs(vel_edge(i, j, k)) * dt / dx;

    if (vel_edge(i, j, k + 1) > small_vel) {
        Ip = sp -
             (0.5 * sigmap) * ((sp - sm) - (1.e0 - 2.e0 / 3.e0 * sigmap) * s6);
    } else {
        Ip = S(i, j, k, n);
    }

    if (vel_edge(i, j, k) < -small_vel) {
        Im = sm +
             (0.5 * sigmam) * ((sp - sm) + (1.e0 - 2.e0 / 3.e0 * sigmam) * s6);
    } else {
        Im = S(i, j, k, n);
    }
}

} // namespace

#endif
