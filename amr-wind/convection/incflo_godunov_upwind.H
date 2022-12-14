#ifndef GODUNOV_UPWIND_H
#define GODUNOV_UPWIND_H

#include <AMReX_Gpu.H>
#include <AMReX_BCRec.H>

/* This header file contains the inlined __host__ __device__ functions required
   for the scalar advection routines for 3D Godunov. It also contains function
   declarations for controlling host functions.  */

namespace {

// This version is called after the MAC projection, when we use the
// MAC-projected velocity for upwinding
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_upwind_fpu(
    const int i,
    const int j,
    const int k,
    const int n,
    amrex::Real& Im,
    amrex::Real& Ip,
    const amrex::Array4<const amrex::Real>& S)
{
    Ip = S(i, j, k, n);
    Im = S(i, j, k, n);
}

} // namespace

#endif
