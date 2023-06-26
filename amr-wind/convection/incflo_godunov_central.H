#ifndef GODUNOV_CENTRAL_H
#define GODUNOV_CENTRAL_H

#include <AMReX_Gpu.H>
#include <AMReX_BCRec.H>

/* This header file contains the inlined __host__ __device__ functions required
   for the scalar advection routines for 3D Godunov. It also contains function
   declarations for controlling host functions.  */

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Godunov_central_fpu(
    const int i,
    const int j,
    const int k,
    const int n,
    amrex::Real& Im,
    amrex::Real& Ip,
    const amrex::Array4<const amrex::Real>& S,
    const int dir)
{
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));
    const amrex::IntVect ivp(iv + amrex::IntVect::TheDimensionVector(dir));
    const amrex::Real sp = S(ivp, n);
    const amrex::Real sm = S(ivm, n);
    Ip = 0.5 * (sp + sm);
    Im = 0.5 * (sp + sm);
}

} // namespace

#endif
