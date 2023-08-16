#ifndef GODUNOV_MINMOD_H
#define GODUNOV_MINMOD_H

#include <AMReX_Gpu.H>
#include <AMReX_BCRec.H>

/* This header file contains the inlined __host__ __device__ functions required
   for the scalar advection routines for 3D Godunov. It also contains function
   declarations for controlling host functions.  */

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void minmod(
    const amrex::Real sm1,
    const amrex::Real s0,
    const amrex::Real sp1,
    amrex::Real& dsm,
    amrex::Real& dsp)
{
    // Calculate gradients on both sides
    dsp = sp1 - s0;
    dsm = s0 - sm1;

    if (!(dsp * dsm < 0.0)) {
        // Select the smaller slope if same sign
        if (std::abs(dsp) < std::abs(dsm)) {
            dsm = dsp;
        } else {
            dsp = dsm;
        }
    } else {
        // Set to zero if opposite sign
        dsp = 0.0;
        dsm = 0.0;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_minmod_bc(
    const int n,
    const amrex::Real sm1,
    const amrex::Real s0,
    const amrex::Real sp1,
    amrex::Real& dsm,
    amrex::Real& dsp,
    const int bclo,
    const int bchi,
    const int domlo,
    const int domhi)
{
    using namespace amrex;

    if (bclo == BCType::ext_dir || bclo == BCType::hoextrap) {
        if (n == domlo) {
            // Ensure that left-side slope is used unchanged
            dsm = s0 - sm1;
        }
    }

    if (bchi == BCType::ext_dir || bchi == BCType::hoextrap) {
        if (n == domhi) {
            // Ensure that the right-side slope is used unchanged
            dsp = sp1 - s0;
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
    amrex::Real dsp = 0.0;
    amrex::Real dsm = 0.0;
    minmod(sm1, s0, sp1, dsm, dsp);
    Godunov_minmod_bc(
        i, sm1, s0, sp1, dsm, dsp, bc.lo(0), bc.hi(0), domlo, domhi);

    amrex::Real sigmap = std::abs(vel_edge(i + 1, j, k)) * dt / dx;
    amrex::Real sigmam = std::abs(vel_edge(i, j, k)) * dt / dx;

    if (vel_edge(i + 1, j, k) > small_vel) {
        Ip = s0 + 0.5 * (1.0 - sigmap) * dsp;
    } else {
        Ip = S(i, j, k, n);
    }

    if (vel_edge(i, j, k) < -small_vel) {
        Im = s0 - 0.5 * (1.0 - sigmam) * dsm;
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
    amrex::Real dsp = 0.0;
    amrex::Real dsm = 0.0;
    minmod(sm1, s0, sp1, dsm, dsp);
    Godunov_minmod_bc(
        j, sm1, s0, sp1, dsm, dsp, bc.lo(1), bc.hi(1), domlo, domhi);

    amrex::Real sigmap = std::abs(vel_edge(i, j + 1, k)) * dt / dx;
    amrex::Real sigmam = std::abs(vel_edge(i, j, k)) * dt / dx;

    if (vel_edge(i, j + 1, k) > small_vel) {
        Ip = s0 + 0.5 * (1.0 - sigmap) * dsp;
    } else {
        Ip = S(i, j, k, n);
    }

    if (vel_edge(i, j, k) < -small_vel) {
        Im = s0 - 0.5 * (1.0 - sigmam) * dsm;
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
    amrex::Real dsp = 0.0;
    amrex::Real dsm = 0.0;
    minmod(sm1, s0, sp1, dsm, dsp);
    Godunov_minmod_bc(
        k, sm1, s0, sp1, dsm, dsp, bc.lo(2), bc.hi(2), domlo, domhi);

    amrex::Real sigmap = std::abs(vel_edge(i, j, k + 1)) * dt / dx;
    amrex::Real sigmam = std::abs(vel_edge(i, j, k)) * dt / dx;

    if (vel_edge(i, j, k + 1) > small_vel) {
        Ip = s0 + 0.5 * (1.0 - sigmap) * dsp;
    } else {
        Ip = S(i, j, k, n);
    }

    if (vel_edge(i, j, k) < -small_vel) {
        Im = s0 - 0.5 * (1.0 - sigmam) * dsm;
    } else {
        Im = S(i, j, k, n);
    }
}

} // namespace

#endif
