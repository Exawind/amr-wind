#ifndef GODUNOV_CENTRAL_H
#define GODUNOV_CENTRAL_H

#include <AMReX_Gpu.H>
#include <AMReX_BCRec.H>

/* This header file contains the inlined __host__ __device__ functions required
   for the scalar advection routines for 3D Godunov. It also contains function
   declarations for controlling host functions.  */

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void central(
    const amrex::Real sm1,
    const amrex::Real s0,
    const amrex::Real sp1,
    amrex::Real& dsm,
    amrex::Real& dsp)
{
    // Calculate gradients on both sides
    dsp = sp1 - s0;
    dsm = s0 - sm1;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_central_bc(
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

// !!! No upwinding for central scheme !!!
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_central_x(
    const int i,
    const int j,
    const int k,
    const int n,
    const amrex::Real dx,
    amrex::Real& Im,
    amrex::Real& Ip,
    const amrex::Array4<const amrex::Real>& S,
    const amrex::BCRec& bc,
    const int domlo,
    const int domhi)
{
    using namespace amrex;
    amrex::Real sm1 = S(i - 1, j, k, n);
    amrex::Real s0 = S(i, j, k, n);
    amrex::Real sp1 = S(i + 1, j, k, n);
    amrex::Real dsp = 0.0;
    amrex::Real dsm = 0.0;
    central(sm1, s0, sp1, dsm, dsp);
    Godunov_central_bc(
        i, sm1, s0, sp1, dsm, dsp, bc.lo(0), bc.hi(0), domlo, domhi);

    // Interpolate to edges
    Ip = s0 + 0.5 * dsp;
    Im = s0 - 0.5 * dsm;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_central_y(
    const int i,
    const int j,
    const int k,
    const int n,
    const amrex::Real dx,
    amrex::Real& Im,
    amrex::Real& Ip,
    const amrex::Array4<const amrex::Real>& S,
    const amrex::BCRec& bc,
    const int domlo,
    const int domhi)
{
    using namespace amrex;
    amrex::Real sm1 = S(i, j - 1, k, n);
    amrex::Real s0 = S(i, j, k, n);
    amrex::Real sp1 = S(i, j + 1, k, n);
    amrex::Real dsp = 0.0;
    amrex::Real dsm = 0.0;
    central(sm1, s0, sp1, dsm, dsp);
    Godunov_central_bc(
        j, sm1, s0, sp1, dsm, dsp, bc.lo(1), bc.hi(1), domlo, domhi);

    Ip = s0 + 0.5 * dsp;
    Im = s0 - 0.5 * dsm;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_central_z(
    const int i,
    const int j,
    const int k,
    const int n,
    const amrex::Real dx,
    amrex::Real& Im,
    amrex::Real& Ip,
    const amrex::Array4<const amrex::Real>& S,
    const amrex::BCRec& bc,
    const int domlo,
    const int domhi)
{
    using namespace amrex;
    amrex::Real sm1 = S(i, j, k - 1, n);
    amrex::Real s0 = S(i, j, k, n);
    amrex::Real sp1 = S(i, j, k + 1, n);
    amrex::Real dsp = 0.0;
    amrex::Real dsm = 0.0;
    central(sm1, s0, sp1, dsm, dsp);
    Godunov_central_bc(
        k, sm1, s0, sp1, dsm, dsp, bc.lo(2), bc.hi(2), domlo, domhi);

    Ip = s0 + 0.5 * dsp;
    Im = s0 - 0.5 * dsm;
}

} // namespace

#endif
