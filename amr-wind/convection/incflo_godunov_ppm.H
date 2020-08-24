#ifndef GODUNOV_PPM_H
#define GODUNOV_PPM_H

#include <AMReX_Gpu.H> 
#include <AMReX_BCRec.H>

/* This header file contains the inlined __host__ __device__ functions required for 
   the scalar advection routines for 3D Godunov. It also contains function declarations 
   for controlling host functions.  */

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real vanLeer(const amrex::Real a, const amrex::Real b, const amrex::Real c)
{
    constexpr amrex::Real small_qty_sq = (1.e-10)*(1.e-10);

    amrex::Real dsc = 0.5e0*(b - c);
    amrex::Real dsl = 2.e0*(a - c);
    amrex::Real dsr = 2.e0*(b - a);
    return (dsl*dsr > small_qty_sq) ? 
        amrex::Math::copysign(1., dsc)*amrex::min(amrex::Math::abs(dsc),amrex::min(amrex::Math::abs(dsl), amrex::Math::abs(dsr))) : 0.;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_ppm_xbc(const int i, const int j, const int k, const int n,
                     amrex::Real &sm, amrex::Real &sp, 
                     amrex::Real &sedge1, amrex::Real &sedge2, 
                     const amrex::Array4<const amrex::Real> &s, 
                     const int bclo, const int bchi, 
                     const int domlo, const int domhi)
{
    using namespace amrex;

    if (bclo == BCType::ext_dir || bclo == BCType::hoextrap)
    {
        if (i == domlo)
        {
            sedge2 = -0.2*s(domlo-1,j,k,n) + 0.75*s(domlo,j, k, n) 
                   +  0.5*s(domlo+1,j,k,n) - 0.05*s(domlo+2,j,k,n);
            sedge2 = amrex::max(sedge2, amrex::min(s(domlo+1,j,k,n), s(domlo,j,k,n)));
            sedge2 = amrex::min(sedge2, amrex::max(s(domlo+1,j,k,n), s(domlo,j,k,n)));

            sm = s(domlo-1,j,k,n);
            sp = sedge2;

        } else if (i == domlo+1) {

            sedge1 = -0.2*s(domlo-1,j,k,n) + 0.75*s(domlo  ,j,k,n) 
                   +  0.5*s(domlo+1,j,k,n) - 0.05*s(domlo+2,j,k,n);
            sedge1 = amrex::max(sedge1, amrex::min(s(domlo+1,j,k,n), s(domlo,j,k,n)));
            sedge1 = amrex::min(sedge1, amrex::max(s(domlo+1,j,k,n), s(domlo,j,k,n)));

            sp = sedge2;
            sm = sedge1;

            if ( (sp - s(domlo+1,j,k,n))*(s(domlo+1,j,k,n) - sm) <= 0.)
            {
                sp = s(domlo+1,j,k,n);
                sm = s(domlo+1,j,k,n);
            }
            else if(amrex::Math::abs(sp - s(domlo+1,j,k,n)) >= 2.*amrex::Math::abs(sm - s(domlo+1,j,k,n)))
                sp = 3.*s(domlo+1,j,k,n) - 2.*sm;
            else if(amrex::Math::abs(sm - s(domlo+1,j,k,n)) >= 2.*amrex::Math::abs(sp - s(domlo+1,j,k,n)))
                sm = 3.*s(domlo+1,j,k,n) - 2.*sp;
        }
    }

    if (bchi == BCType::ext_dir || bchi == BCType::hoextrap)
    {
        if (i == domhi)
        { 
            sedge1 =  -0.2*s(domhi+1,j,k,n) + 0.75*s(domhi,j,k, n) 
                      +0.5*s(domhi-1,j,k,n) - 0.05*s(domhi-2,j,k,n);
            sedge1 = amrex::max(sedge1, amrex::min(s(domhi-1,j,k,n), s(domhi,j,k,n)));
            sedge1 = amrex::min(sedge1, amrex::max(s(domhi-1,j,k,n), s(domhi,j,k,n)));

            sp = s(domhi+1,j,k,n);
            sm = sedge1;

        } else if (i == domhi-1) {

            sedge2 = -0.2*s(domhi+1,j,k,n) + 0.75*s(domhi  ,j,k,n) 
                     +0.5*s(domhi-1,j,k,n) - 0.05*s(domhi-2,j,k,n);
            sedge2 = amrex::max(sedge2, amrex::min(s(domhi-1,j,k,n), s(domhi,j,k,n)));
            sedge2 = amrex::min(sedge2, amrex::max(s(domhi-1,j,k,n), s(domhi,j,k,n)));

            sp = sedge2;
            sm = sedge1;

            if( (sp - s(domhi-1,j,k,n))*(s(domhi-1,j,k,n) - sm) <= 0.)
            {
                sp = s(domhi-1,j,k,n);
                sm = s(domhi-1,j,k,n);
            }
            else if(amrex::Math::abs(sp - s(domhi-1,j,k,n)) >= 2.*amrex::Math::abs(sm - s(domhi-1,j,k,n)))
                sp = 3.*s(domhi-1,j,k,n) - 2.*sm;
            else if(amrex::Math::abs(sm - s(domhi-1,j,k,n)) >= 2.*amrex::Math::abs(sp - s(domhi-1,j,k,n)))
                sm = 3.*s(domhi-1,j,k,n) - 2.*sp;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_ppm_ybc(const int i, const int j, const int k, const int n,
                           amrex::Real &sm, amrex::Real &sp, 
                           amrex::Real &sedge1, amrex::Real &sedge2, 
                           const amrex::Array4<const amrex::Real> &s, 
                           const int bclo, const int bchi, 
                           const int domlo, const int domhi)
{
    using namespace amrex;

    if (bclo == BCType::ext_dir || bclo == BCType::hoextrap)
    {
        if (j == domlo)
        {
            sedge2 = -0.2*s(i,domlo-1,k,n) + 0.75*s(i,domlo  ,k,n) 
                     +0.5*s(i,domlo+1,k,n) - 0.05*s(i,domlo+2,k,n);
            sedge2 = amrex::max(sedge2, amrex::min(s(i,domlo+1,k,n), s(i,domlo,k,n)));
            sedge2 = amrex::min(sedge2, amrex::max(s(i,domlo+1,k,n), s(i,domlo,k,n)));

            sm = s(i,domlo-1,k,n);
            sp = sedge2;

        } else if (j == domlo+1) {

            sedge1 = -0.2*s(i,domlo-1,k,n) + 0.75*s(i,domlo  ,k,n) 
                     +0.5*s(i,domlo+1,k,n) - 0.05*s(i,domlo+2,k,n);
            sedge1 = amrex::max(sedge1, amrex::min(s(i,domlo+1,k,n), s(i,domlo,k,n)));
            sedge1 = amrex::min(sedge1, amrex::max(s(i,domlo+1,k,n), s(i,domlo,k,n)));

            sp = sedge2;
            sm = sedge1;

            if ( (sp - s(i,domlo+1,k,n))*(s(i,domlo+1,k,n) - sm) <= 0.)
            {
                sp = s(i,domlo+1,k,n);
                sm = s(i,domlo+1,k,n);
            }
            else if(amrex::Math::abs(sp - s(i,domlo+1,k,n)) >= 2.*amrex::Math::abs(sm - s(i,domlo+1,k,n)))
                sp = 3.*s(i,domlo+1,k,n) - 2.*sm;
            else if(amrex::Math::abs(sm - s(i,domlo+1,k,n)) >= 2.*amrex::Math::abs(sp - s(i,domlo+1,k,n)))
                sm = 3.*s(i,domlo+1,k,n) - 2.*sp;
        }
    }

    if (bchi == BCType::ext_dir || bchi == BCType::hoextrap)
    {
        if (j == domhi)
        {
            sedge1 =  -0.2*s(i,domhi+1,k,n) + 0.75*s(i,domhi  ,k,n) 
                      +0.5*s(i,domhi-1,k,n) - 0.05*s(i,domhi-2,k,n);
            sedge1 = amrex::max(sedge1, amrex::min(s(i,domhi-1,k,n), s(i,domhi,k,n)));
            sedge1 = amrex::min(sedge1, amrex::max(s(i,domhi-1,k,n), s(i,domhi,k,n)));

            sp = s(i,domhi+1,k, n);
            sm = sedge1;

        } else if (j == domhi-1) {

            sedge2 =  -0.2*s(i,domhi+1,k,n) + 0.75*s(i,domhi  ,k,n) 
                      +0.5*s(i,domhi-1,k,n) - 0.05*s(i,domhi-2,k,n);
            sedge2 = amrex::max(sedge2, amrex::min(s(i,domhi-1,k,n), s(i,domhi,k,n)));
            sedge2 = amrex::min(sedge2, amrex::max(s(i,domhi-1,k,n), s(i,domhi,k,n)));

            sp = sedge2;
            sm = sedge1;

            if( (sp - s(i,domhi-1,k,n))*(s(i,domhi-1,k,n) - sm) <= 0.){
                sp = s(i,domhi-1,k,n);
                sm = s(i,domhi-1,k,n);
            }
            else if(amrex::Math::abs(sp - s(i,domhi-1,k,n)) >= 2.*amrex::Math::abs(sm - s(i,domhi-1,k,n))) 
                sp = 3.*s(i,domhi-1,k,n) - 2.*sm;

            else if(amrex::Math::abs(sm - s(i,domhi-1,k,n)) >= 2.*amrex::Math::abs(sp - s(i,domhi-1,k,n))) 
                sm = 3.*s(i,domhi-1,k,n) - 2.*sp;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_ppm_zbc(const int i, const int j, const int k, const int n,
                           amrex::Real &sm, amrex::Real &sp, 
                           amrex::Real &sedge1, amrex::Real &sedge2, 
                           const amrex::Array4<const amrex::Real> &s, 
                           const int bclo, const int bchi, 
                           const int domlo, const int domhi)
{
    using namespace amrex;

    if (bclo == BCType::ext_dir || bclo == BCType::hoextrap)
    {

        if (k == domlo)
        {
            sedge2 = -0.2*s(i,j,domlo-1,n) + 0.75*s(i,j,domlo  ,n) 
                     +0.5*s(i,j,domlo+1,n) - 0.05*s(i,j,domlo+2,n);
            sedge2 = amrex::max(sedge2, amrex::min(s(i,j,domlo+1,n), s(i,j,domlo,n)));
            sedge2 = amrex::min(sedge2, amrex::max(s(i,j,domlo+1,n), s(i,j,domlo,n)));

            sm = s(i,j,domlo-1,n);
            sp = sedge2;

        } else if (k == domlo+1) {

            sedge1 = -0.2*s(i,j,domlo-1,n) + 0.75*s(i,j,domlo  ,n) 
                     +0.5*s(i,j,domlo+1,n) - 0.05*s(i,j,domlo+2,n);
            sedge1 = amrex::max(sedge1, amrex::min(s(i,j,domlo+1,n), s(i,j,domlo,n)));
            sedge1 = amrex::min(sedge1, amrex::max(s(i,j,domlo+1,n), s(i,j,domlo,n)));

            sp = sedge2;
            sm = sedge1;

            if ( (sp - s(i,j,domlo+1,n))*(s(i,j,domlo+1,n) - sm) <= 0. )
            {
                sp = s(i,j,domlo+1,n);
                sm = s(i,j,domlo+1,n);
            }
            else if(amrex::Math::abs(sp - s(i,j,domlo+1,n)) >= 2.*amrex::Math::abs(sm - s(i,j,domlo+1,n))) 
                sp = 3.*s(i,j,domlo+1,n) - 2.*sm;

            else if(amrex::Math::abs(sm - s(i,j,domlo+1,n)) >= 2.*amrex::Math::abs(sp - s(i,j,domlo+1,n)))
                sm = 3.*s(i,j,domlo+1,n) - 2.*sp;
        }
    }

    if (bchi == BCType::ext_dir || bchi == BCType::hoextrap)
    {
        if (k == domhi) 
        { 
            sedge1 =  -0.2*s(i,j,domhi+1,n) + 0.75*s(i,j,domhi  ,n) 
                      +0.5*s(i,j,domhi-1,n) - 0.05*s(i,j,domhi-2,n);
            sedge1 = amrex::max(sedge1, amrex::min(s(i,j,domhi-1,n), s(i,j,domhi,n)));
            sedge1 = amrex::min(sedge1, amrex::max(s(i,j,domhi-1,n), s(i,j,domhi,n)));

            sp = s(i,j,domhi+1,n);
            sm = sedge1;

        } else if (k == domhi-1) {

            sedge2 = -0.2*s(i,j,domhi+1,n) + 0.75*s(i,j,domhi  ,n) 
                     +0.5*s(i,j,domhi-1,n) - 0.05*s(i,j,domhi-2,n);
            sedge2 = amrex::max(sedge2, amrex::min(s(i,j,domhi-1,n), s(i,j,domhi,n)));
            sedge2 = amrex::min(sedge2, amrex::max(s(i,j,domhi-1,n), s(i,j,domhi,n)));

            sp = sedge2;
            sm = sedge1;

            if ( (sp - s(i,j,domhi-1,n))*(s(i,j,domhi-1,n) - sm) <= 0. )
            {
                sp = s(i,j,domhi-1,n);
                sm = s(i,j,domhi-1,n);
            }
            else if(amrex::Math::abs(sp - s(i,j,domhi-1,n)) >= 2.*amrex::Math::abs(sm - s(i,j,domhi-1,n)))
                sp = 3.*s(i,j,domhi-1,n) - 2.*sm;
            
            else if(amrex::Math::abs(sm - s(i,j,domhi-1,n)) >= 2.*amrex::Math::abs(sp - s(i,j,domhi-1,n))) 
                sm = 3.*s(i,j,domhi-1,n) - 2.*sp;
              
        }
    }
} 

// Right now only ppm type 1 is supported on GPU
// This version is called before the MAC projection, when we use the cell-centered velocity
//      for upwinding
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_ppm_pred_x(const int i, const int j, const int k, const int n,
                        const amrex::Real dtdx, const amrex::Real v_ad, 
                        const amrex::Array4<const amrex::Real> &S, 
                        const amrex::Array4<amrex::Real> &Im, 
                        const amrex::Array4<amrex::Real> &Ip,
                        const amrex::BCRec bc, const int domlo, const int domhi)
{
    using namespace amrex;

    constexpr amrex::Real small_vel = 1e-10;
    constexpr amrex::Real sixth = 1.0/6.0;

    amrex::Real sedge1, sedge2, sp, sm;
 
    amrex::Real sm2 = S(i-2,j,k,n);
    amrex::Real sm1 = S(i-1,j,k,n);
    amrex::Real s0  = S(i  ,j,k,n);
    amrex::Real sp1 = S(i+1,j,k,n);
    amrex::Real sp2 = S(i+2,j,k,n);

    amrex::Real d1 = vanLeer(s0,sp1,sm1);
    amrex::Real d2 = vanLeer(sm1,s0,sm2);

    sedge1 = 0.5e0*(s0 + sm1) - sixth*(d1 - d2);
    sedge1 = amrex::min(amrex::max(sedge1, amrex::min(s0, sm1)),amrex::max(s0,sm1));

    d1 = vanLeer(sp1,sp2,s0);
    d2 = vanLeer(s0,sp1,sm1);

    sedge2 = 0.5e0*(sp1 + s0) - sixth*(d1 - d2);
    sedge2 = amrex::min(amrex::max(sedge2, amrex::min(s0, sp1)),amrex::max(s0,sp1));
    
    sm = sedge1;
    sp = sedge2;

    if ((sedge2-s0)*(s0-sedge1) < 0.e0)
    {
        sp = s0;
        sm = s0;

    } else if (amrex::Math::abs(sedge2-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge1-s0))
        sp = 3.0*s0 - 2.0*sedge1;

      else if (amrex::Math::abs(sedge1-S(i,j,k,n)) >=  2.0*amrex::Math::abs(sedge2-s0))
        sm = 3.0*s0 - 2.0*sedge2;

    Godunov_ppm_xbc(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(0), bc.hi(0), domlo, domhi);

    amrex::Real s6 = 6.0*s0 - 3.0*(sm + sp);

    amrex::Real sigma = amrex::Math::abs(v_ad)*dtdx;

    if (v_ad > small_vel)
    {
        Ip(i,j,k,n) = sp - (0.5*sigma)*((sp-sm) - (1.0 - 2.0/3.0*sigma)*s6);
        Im(i,j,k,n) = S(i,j,k,n);
    }
    else if (v_ad < -small_vel)
    {
        Ip(i,j,k,n) = S(i,j,k,n);
        Im(i,j,k,n) = sm + (0.5*sigma)*((sp-sm) + (1.0 - 2.0/3.0*sigma)*s6);
    } else 
    {
        Ip(i,j,k,n) = S(i,j,k,n);
        Im(i,j,k,n) = S(i,j,k,n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_ppm_pred_y(const int i, const int j, const int k, const int n,
                        const amrex::Real dtdy, const amrex::Real v_ad, 
                        const amrex::Array4<const amrex::Real> &S,
                        const amrex::Array4<amrex::Real> &Im,
                        const amrex::Array4<amrex::Real> &Ip,
                        const amrex::BCRec bc, const int domlo, const int domhi)
{
    using namespace amrex;

    constexpr amrex::Real small_vel = 1e-10;
    constexpr amrex::Real sixth = 1.0/6.0;

    amrex::Real sedge1, sedge2, sp, sm;
 
    amrex::Real sm2 = S(i,j-2,k,n);
    amrex::Real sm1 = S(i,j-1,k,n);
    amrex::Real s0  = S(i,j  ,k,n);
    amrex::Real sp1 = S(i,j+1,k,n);
    amrex::Real sp2 = S(i,j+2,k,n);

    amrex::Real d1 = vanLeer(s0,sp1,sm1);
    amrex::Real d2 = vanLeer(sm1,s0,sm2);

    sedge1 = 0.5e0*(s0 + sm1) - sixth*(d1 - d2);
    sedge1 = amrex::min(amrex::max(sedge1, amrex::min(s0, sm1)),amrex::max(s0,sm1));

    d1 = vanLeer(sp1,sp2,s0);
    d2 = vanLeer(s0,sp1,sm1);

    sedge2 = 0.5e0*(sp1 + s0) - sixth*(d1 - d2);
    sedge2 = amrex::min(amrex::max(sedge2, amrex::min(s0, sp1)),amrex::max(s0,sp1));
    
    sm = sedge1;
    sp = sedge2;

    if ((sedge2-s0)*(s0-sedge1) < 0.e0)
    {
        sp = s0;
        sm = s0;

    } else if (amrex::Math::abs(sedge2-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge1-s0))
        sp = 3.0*s0 - 2.0*sedge1;

      else if (amrex::Math::abs(sedge1-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge2-s0))
        sm = 3.0*s0 - 2.0*sedge2;

    Godunov_ppm_ybc(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(1), bc.hi(1), domlo, domhi);

    amrex::Real s6 = 6.0*s0- 3.0*(sm + sp);

    amrex::Real sigma = amrex::Math::abs(v_ad)*dtdy;

    if (v_ad > small_vel)
    {
        Ip(i,j,k,n) = sp - (0.5*sigma)*((sp-sm) - (1.0 - 2.0/3.0*sigma)*s6);
        Im(i,j,k,n) = S(i,j,k,n);
    } 
    else if (v_ad < -small_vel)
    {
        Ip(i,j,k,n) = S(i,j,k,n);
        Im(i,j,k,n) = sm + (0.5*sigma)*((sp-sm) + (1.0 - 2.0/3.0*sigma)*s6);
    } else 
    {
        Ip(i,j,k,n) = S(i,j,k,n);
        Im(i,j,k,n) = S(i,j,k,n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_ppm_pred_z(const int i, const int j, const int k, const int n,
                        const amrex::Real dtdz, const amrex::Real v_ad, 
                        const amrex::Array4<const amrex::Real> &S,
                        const amrex::Array4<amrex::Real> &Im,
                        const amrex::Array4<amrex::Real> &Ip,
                        const amrex::BCRec bc, const int domlo, const int domhi)
{
    using namespace amrex;

    constexpr amrex::Real small_vel = 1e-10;
    constexpr amrex::Real sixth = 1.0/6.0;

    amrex::Real sedge1, sedge2, sp, sm;
 
    amrex::Real sm2 = S(i,j,k-2,n);
    amrex::Real sm1 = S(i,j,k-1,n);
    amrex::Real s0  = S(i,j,k  ,n);
    amrex::Real sp1 = S(i,j,k+1,n);
    amrex::Real sp2 = S(i,j,k+2,n);

    amrex::Real d1 = vanLeer(s0,sp1,sm1);
    amrex::Real d2 = vanLeer(sm1,s0,sm2);

    sedge1 = 0.5e0*(s0 + sm1) - sixth*(d1 - d2);
    sedge1 = amrex::min(amrex::max(sedge1, amrex::min(s0, sm1)),amrex::max(s0,sm1));

    d1 = vanLeer(sp1,sp2,s0);
    d2 = vanLeer(s0,sp1,sm1);

    sedge2 = 0.5e0*(sp1 + s0) - sixth*(d1 - d2);
    sedge2 = amrex::min(amrex::max(sedge2, amrex::min(s0, sp1)),amrex::max(s0,sp1));
    
    sm = sedge1;
    sp = sedge2;

    if ((sedge2-s0)*(s0-sedge1) < 0.e0)
    {
        sp = s0;
        sm = s0;

    } else if (amrex::Math::abs(sedge2-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge1-s0))
        sp = 3.0*s0 - 2.0*sedge1;

      else if (amrex::Math::abs(sedge1-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge2-s0))
        sm = 3.0*s0 - 2.0*sedge2;

    Godunov_ppm_zbc(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(2), bc.hi(2), domlo, domhi);

    amrex::Real s6 = 6.0*s0- 3.0*(sm + sp);

    amrex::Real sigma = amrex::Math::abs(v_ad)*dtdz;

    if (v_ad > small_vel)
    {
        Ip(i,j,k,n) = sp - (0.5*sigma)*((sp-sm) - (1.0 - 2.0/3.0*sigma)*s6);
        Im(i,j,k,n) = S(i,j,k,n);
    }
    else if (v_ad < -small_vel)
    {
        Ip(i,j,k,n) = S(i,j,k,n);
        Im(i,j,k,n) = sm + (0.5*sigma)*((sp-sm) + (1.0 - 2.0/3.0*sigma)*s6);
    } else 
    {
        Ip(i,j,k,n) = S(i,j,k,n);
        Im(i,j,k,n) = S(i,j,k,n);
    }
}

// Right now only ppm type 1 is supported on GPU
// This version is called after the MAC projection, when we use the MAC-projected velocity
//      for upwinding
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_ppm_fpu_x (const int i, const int j, const int k, const int n,
                        const amrex::Real dt, const amrex::Real dx,
                        amrex::Real& Im, amrex::Real& Ip,
                        const amrex::Array4<const amrex::Real> &S,
                        const amrex::Array4<const amrex::Real> &vel_edge,
                        const amrex::BCRec bc,
                        const int domlo, const int domhi)
{
    {
    using namespace amrex;

    constexpr Real small_vel = 1e-10;
    constexpr Real sixth = 1.0/6.0;

    amrex::Real sedge1, sedge2, sp, sm;
 
    amrex::Real sm2 = S(i-2,j,k,n);
    amrex::Real sm1 = S(i-1,j,k,n);
    amrex::Real s0  = S(i  ,j,k,n);
    amrex::Real sp1 = S(i+1,j,k,n);
    amrex::Real sp2 = S(i+2,j,k,n);

    amrex::Real d1 = vanLeer(s0,sp1,sm1);
    amrex::Real d2 = vanLeer(sm1,s0,sm2);

    sedge1 = 0.5e0*(s0 + sm1) - sixth*(d1 - d2);
    sedge1 = amrex::min(amrex::max(sedge1, amrex::min(s0, sm1)),amrex::max(s0,sm1));

    d1 = vanLeer(sp1,sp2,s0);
    d2 = vanLeer(s0,sp1,sm1);

    sedge2 = 0.5e0*(sp1 + s0) - sixth*(d1 - d2);
    sedge2 = amrex::min(amrex::max(sedge2, amrex::min(s0, sp1)),amrex::max(s0,sp1));
    
    sm = sedge1;
    sp = sedge2;

    if ((sedge2-s0)*(s0-sedge1) < 0.e0)
    {
        sp = s0;
        sm = s0;
    }
    else if (amrex::Math::abs(sedge2-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge1-s0))
      sp = 3.0*s0 - 2.0*sedge1;

    else if (amrex::Math::abs(sedge1-S(i,j,k,n)) >=  2.0*amrex::Math::abs(sedge2-s0))
      sm = 3.0*s0 - 2.0*sedge2;

    Godunov_ppm_xbc(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(0), bc.hi(0), domlo, domhi);

    Real s6 = 6.0*s0 - 3.0*(sm + sp);

    Real sigmap = amrex::Math::abs(vel_edge(i+1,j,k))*dt/dx;
    Real sigmam = amrex::Math::abs(vel_edge(i  ,j,k))*dt/dx;

    if (vel_edge(i+1,j,k) > small_vel) 
        Ip = sp - (0.5*sigmap)*((sp - sm) - (1.e0 -2.e0/3.e0*sigmap)*s6);
    else
        Ip = S(i,j,k,n);

    if(vel_edge(i,j,k) < -small_vel)
        Im = sm + (0.5*sigmam)*((sp-sm) + (1.e0 - 2.e0/3.e0*sigmam)*s6);
    else
        Im = S(i,j,k,n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_ppm_fpu_y (const int i, const int j, const int k, const int n,
                        const amrex::Real dt, const amrex::Real dx,
                        amrex::Real& Im, amrex::Real& Ip,
                        const amrex::Array4<const amrex::Real> &S,
                        const amrex::Array4<const amrex::Real> &vel_edge,
                        const amrex::BCRec bc,
                        const int domlo, const int domhi)
{
    {
    using namespace amrex;

    constexpr amrex::Real small_vel = 1e-10;
    constexpr amrex::Real sixth = 1.0/6.0;

    Real sedge1, sedge2, sp, sm;
 
    Real sm2 = S(i,j-2,k,n);
    Real sm1 = S(i,j-1,k,n);
    Real s0  = S(i,j  ,k,n);
    Real sp1 = S(i,j+1,k,n);
    Real sp2 = S(i,j+2,k,n);

    Real d1 = vanLeer(s0,sp1,sm1);
    Real d2 = vanLeer(sm1,s0,sm2);

    sedge1 = 0.5e0*(s0 + sm1) - sixth*(d1 - d2);
    sedge1 = amrex::min(amrex::max(sedge1, amrex::min(s0, sm1)),amrex::max(s0,sm1));

    d1 = vanLeer(sp1,sp2,s0);
    d2 = vanLeer(s0,sp1,sm1);

    sedge2 = 0.5e0*(sp1 + s0) - sixth*(d1 - d2);
    sedge2 = amrex::min(amrex::max(sedge2, amrex::min(s0, sp1)),amrex::max(s0,sp1));
    
    sm = sedge1;
    sp = sedge2;

    if ((sedge2-s0)*(s0-sedge1) < 0.e0)
    {
        sp = s0;
        sm = s0;
    } 
    else if (amrex::Math::abs(sedge2-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge1-s0))
        sp = 3.0*s0 - 2.0*sedge1;

    else if (amrex::Math::abs(sedge1-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge2-s0))
        sm = 3.0*s0 - 2.0*sedge2;

    Godunov_ppm_ybc(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(1), bc.hi(1), domlo, domhi);

    amrex::Real s6 = 6.0*s0- 3.0*(sm + sp);

    Real sigmap = amrex::Math::abs(vel_edge(i,j+1,k))*dt/dx;
    Real sigmam = amrex::Math::abs(vel_edge(i,j  ,k))*dt/dx;

    if (vel_edge(i,j+1,k) > small_vel)
        Ip = sp - (0.5*sigmap)*((sp - sm) - (1.e0 -2.e0/3.e0*sigmap)*s6);
    else
        Ip = S(i,j,k,n);

    if (vel_edge(i,j,k) < -small_vel)
        Im = sm + (0.5*sigmam)*((sp-sm) + (1.e0 - 2.e0/3.e0*sigmam)*s6);
    else
        Im = S(i,j,k,n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_ppm_fpu_z (const int i, const int j, const int k, const int n,
                        const amrex::Real dt, const amrex::Real dx,
                        amrex::Real& Im, amrex::Real& Ip,
                        const amrex::Array4<const amrex::Real> &S,
                        const amrex::Array4<const amrex::Real> &vel_edge,
                        const amrex::BCRec bc,
                        const int domlo, const int domhi)
{
    {
    using namespace amrex;

    constexpr amrex::Real small_vel = 1e-10;
    constexpr amrex::Real sixth = 1.0/6.0;

    Real sedge1, sedge2, sp, sm;
 
    Real sm2 = S(i,j,k-2,n);
    Real sm1 = S(i,j,k-1,n);
    Real s0  = S(i,j,k  ,n);
    Real sp1 = S(i,j,k+1,n);
    Real sp2 = S(i,j,k+2,n);

    Real d1 = vanLeer(s0,sp1,sm1);
    Real d2 = vanLeer(sm1,s0,sm2);

    sedge1 = 0.5e0*(s0 + sm1) - sixth*(d1 - d2);
    sedge1 = amrex::min(amrex::max(sedge1, amrex::min(s0, sm1)),amrex::max(s0,sm1));

    d1 = vanLeer(sp1,sp2,s0);
    d2 = vanLeer(s0,sp1,sm1);

    sedge2 = 0.5e0*(sp1 + s0) - sixth*(d1 - d2);
    sedge2 = amrex::min(amrex::max(sedge2, amrex::min(s0, sp1)),amrex::max(s0,sp1));
    
    sm = sedge1;
    sp = sedge2;

    if ((sedge2-s0)*(s0-sedge1) < 0.e0)
    {
        sp = s0;
        sm = s0;
    } 
    else if (amrex::Math::abs(sedge2-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge1-s0))
        sp = 3.0*s0 - 2.0*sedge1;

    else if (amrex::Math::abs(sedge1-S(i,j,k,n)) >= 2.0*amrex::Math::abs(sedge2-s0))
        sm = 3.0*s0 - 2.0*sedge2;

    Godunov_ppm_zbc(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(2), bc.hi(2), domlo, domhi);

    Real s6 = 6.0*s0- 3.0*(sm + sp);
    Real sigmap = amrex::Math::abs(vel_edge(i,j,k+1))*dt/dx;
    Real sigmam = amrex::Math::abs(vel_edge(i,j,k  ))*dt/dx;

    if(vel_edge(i,j,k+1) > small_vel)
        Ip = sp - (0.5*sigmap)*((sp-sm) - (1.e0 -2.e0/3.e0*sigmap)*s6);
    else
        Ip = S(i,j,k,n);

    if(vel_edge(i,j,k) < -small_vel)
        Im = sm + (0.5*sigmam)*((sp-sm) + (1.e0 - 2.e0/3.e0*sigmam)*s6);
    else
        Im = S(i,j,k,n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_trans_xbc(const int i, const int j, const int k, const int n,
                              const amrex::Array4<const amrex::Real> &s, 
                                    amrex::Real &lo,
                                    amrex::Real &hi,
                                    amrex::Real & /* uedge */,
                              const int bclo, const int bchi, 
                              const int domlo, const int domhi)
{
    using namespace amrex;

    // constexpr amrex::Real small_vel = 1e-10;

    // Low X
    if (i <= domlo)
    {
        if (bclo==BCType::ext_dir)
        {
            // IAMR does this but it breaks lo/hi symmetry
            // Real st = (uedge <= small_vel) ? hi : s(domlo-1,j,k,n);
            // So here we do something simpler...
            Real st = s(domlo-1,j,k,n);
            lo = st;
            hi = st;
        } 

        else if(bclo == BCType::foextrap || bclo == BCType::hoextrap || bclo == BCType::reflect_even)
            lo = hi;

        else if(bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }

    // High X
    if (i > domhi) 
    {
        if (bchi==BCType::ext_dir)
        {
            // IAMR does this but it breaks lo/hi symmetry
            // Real st = (uedge >= -small_vel)? lo : s(domhi+1,j,k,n);
            // So here we do something simpler...
            Real st = s(domhi+1,j,k,n);
            lo = st;
            hi = st;
        }

        else if (bchi == BCType::foextrap || bchi == BCType::hoextrap || bchi == BCType::reflect_even)
            hi = lo;

        else if(bchi == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_trans_ybc(const int i, const int j, const int k, const int n,
                       const amrex::Array4<const amrex::Real> &s,
                       amrex::Real &lo,
                       amrex::Real &hi,
                       amrex::Real /* vedge */,
                       const int bclo, const int bchi,
                       const int domlo, const int domhi)
{
    using namespace amrex;

    // constexpr amrex::Real small_vel = 1e-10;

    // Low Y 
    if (j <= domlo)
    {
        if (bclo==BCType::ext_dir)
        {
            // IAMR does this but it breaks lo/hi symmetry
            // Real st = (vedge <= small_vel) ? hi : s(i,domlo-1,k,n);
            // So here we do something simpler...
            Real st = s(i,domlo-1,k,n);
            lo = st;
            hi = st;
        }

        else if (bclo == BCType::foextrap || bclo == BCType::hoextrap || bclo == BCType::reflect_even)
            lo = hi;
        
        else if(bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }

    // High Y
    if (j > domhi)
    {
        if (bchi==BCType::ext_dir)
        {
            // IAMR does this but it breaks lo/hi symmetry
            // Real st = (vedge >= -small_vel)? lo : s(i,domhi+1,k,n);
            // So here we do something simpler...
            Real st = s(i,domhi+1,k,n);
            lo = st;
            hi = st;
        }

        else if (bchi == BCType::foextrap || bchi == BCType::hoextrap || bchi == BCType::reflect_even)
            hi = lo;

        else if (bchi == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_trans_zbc(const int i, const int j, const int k, const int n,
                       const amrex::Array4<const amrex::Real> &s,
                       amrex::Real &lo,
                       amrex::Real &hi,
                       amrex::Real /* wedge */,
                       const int bclo, const int bchi,
                       const int domlo, const int domhi)
{
    using namespace amrex;

    // constexpr amrex::Real small_vel = 1e-10;

    // Low Z
    if (k <= domlo)
    { 
        if (bclo==BCType::ext_dir)
        {
            // IAMR does this but it breaks lo/hi symmetry
            // Real st = (wedge <= small_vel) ? hi : s(i,j,domlo-1,n);
            // So here we do something simpler...
            Real st = s(i,j,domlo-1,n);
            lo = st;
            hi = st;
        }

        else if (bclo == BCType::foextrap || bclo == BCType::hoextrap || bclo == BCType::reflect_even)
            lo = hi;

        else if(bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }

    // High Z 
    if (k > domhi)
    {
        if (bchi==BCType::ext_dir)
        {
            // IAMR does this but it breaks lo/hi symmetry
            // Real st = (wedge >= -small_vel)? lo : s(i,j,domhi+1,n);
            // So here we do something simpler...
            Real st = s(i,j,domhi+1,n);
            lo = st;
            hi = st;
        }
        else if (bchi == BCType::foextrap || bchi == BCType::hoextrap || bchi == BCType::reflect_even)
            hi = lo;
        
        else if (bchi == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_cc_xbc_lo(const int i, const int j, const int k, const int n,
                       const amrex::Array4<const amrex::Real> &s, 
                       amrex::Real &lo, amrex::Real &hi, 
                       const amrex::Array4< const amrex::Real> &uedge, 
                       const int bclo, const int domlo)
{
    using namespace amrex;
 
    if (i == domlo)
    {
            if (bclo == BCType::ext_dir && uedge(i, j, k) >= 0.)
            {
                hi = s(domlo-1, j, k, n);
                lo = hi;
            }    
            else if ( bclo == BCType::foextrap || bclo == BCType::hoextrap ||
                      bclo == BCType::reflect_even ||
                     (bclo == BCType::ext_dir && uedge(i, j, k) < 0) )
                lo = hi;

            else if (bclo == BCType::reflect_odd)
            {
                hi = 0.;
                lo = hi;
            }
    } else return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_cc_xbc_hi(const int i, const int j, const int k, const int n,
                       const amrex::Array4<const amrex::Real> &s, 
                       amrex::Real &lo, amrex::Real &hi, 
                       const amrex::Array4< const amrex::Real> &uedge, 
                       const int bchi, const int domhi)
{
    using namespace amrex;

    if (i == domhi+1)
    {
            if (bchi == BCType::ext_dir && uedge(i, j, k) <= 0.)
            {
                lo = s(domhi+1, j, k, n); 
                hi = lo;
            }    
            else if ( bchi == BCType::foextrap || bchi == BCType::hoextrap ||
                      bchi == BCType::reflect_even ||
                     (bchi == BCType::ext_dir && uedge(i, j, k) > 0) )
                hi = lo;

            else if (bchi == BCType::reflect_odd)
            {
                lo = 0.;
                hi = lo;
            }
    }
    else return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_cc_ybc_lo(const int i, const int j, const int k, const int n,
                       const amrex::Array4<const amrex::Real> &s, 
                       amrex::Real &lo, amrex::Real &hi, 
                       const amrex::Array4< const amrex::Real> &vedge, 
                       const int bclo, const int domlo)
{
    using namespace amrex;

    if (j == domlo)
    {
            if (bclo == BCType::ext_dir && vedge(i, j, k) >= 0.)
            {
                hi = s(i, domlo-1, k, n);
                lo = hi;
            }    
            else if ( bclo == BCType::foextrap || bclo == BCType::hoextrap ||
                      bclo == BCType::reflect_even ||
                     (bclo == BCType::ext_dir && vedge(i, j, k) < 0) )
                lo = hi;

            else if (bclo == BCType::reflect_odd){
                hi = 0.;
                lo = hi;
            }
    } else return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_cc_ybc_hi(const int i, const int j, const int k, const int n,
                       const amrex::Array4<const amrex::Real> &s, 
                       amrex::Real &lo, amrex::Real &hi, 
                       const amrex::Array4<const amrex::Real> &vedge, 
                       const int bchi, const int domhi)
{
    using namespace amrex;

    if (j == domhi+1)
    {
            if (bchi == BCType::ext_dir && vedge(i, j, k) <= 0.)
            {
                lo = s(i, domhi+1, k, n);
                hi = lo;
            }    
            else if ( bchi == BCType::foextrap || bchi == BCType::hoextrap ||
                      bchi == BCType::reflect_even ||
                     (bchi == BCType::ext_dir && vedge(i, j, k) > 0) )
                hi = lo;
            
            else if(bchi == BCType::reflect_odd)
            {
                lo = 0.;
                hi = lo;
            }
    }
    else return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_cc_zbc_lo(const int i, const int j, const int k, const int n,
                       const amrex::Array4<const amrex::Real> &s, 
                       amrex::Real &lo, amrex::Real &hi, 
                       const amrex::Array4< const amrex::Real> &wedge, 
                       const int bclo, const int domlo)
{
    using namespace amrex;

    if (k == domlo)
    {
            if (bclo == BCType::ext_dir && wedge(i, j, k) >= 0.)
            {
                hi = s(i, j, domlo-1, n);
                lo = hi;
            }    

            else if ( bclo == BCType::foextrap || bclo == BCType::hoextrap ||
                      bclo == BCType::reflect_even ||
                     (bclo == BCType::ext_dir && wedge(i, j, k) < 0) )
                lo = hi;

            else if(bclo == BCType::reflect_odd)
            {
                hi = 0.;
                lo = hi;
            }
    } else return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_cc_zbc_hi(const int i, const int j, const int k, const int n,
                       const amrex::Array4<const amrex::Real> &s, 
                       amrex::Real &lo, amrex::Real &hi, 
                       const amrex::Array4< const amrex::Real> &wedge, 
                       const int bchi, const int domhi)
{
    using namespace amrex;

    if (k == domhi+1) 
    {
            if (bchi == BCType::ext_dir && wedge(i, j, k) <= 0.)
            {
                lo = s(i,j,domhi+1, n);
                hi = lo;
            }    
            else if ( bchi == BCType::foextrap || bchi == BCType::hoextrap ||
                      bchi == BCType::reflect_even ||
                     (bchi == BCType::ext_dir && wedge(i, j, k) > 0) )
                hi = lo;

            else if(bchi == BCType::reflect_odd)
            {
                lo = 0.;
                hi = lo;
            }
    } else return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_corner_couple_yx (amrex::Real& lo1, amrex::Real& hi1,
                               int i, int j, int k, int n, amrex::Real dt, amrex::Real dx,
                               bool iconserv, amrex::Real lo, amrex::Real hi,
                               amrex::Array4<amrex::Real const> const& s,
                               amrex::Array4<amrex::Real const> const& mac,
                               amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo - dt/(3.e0*dx)*(state(i+1,j-1,k,n)*mac(i+1,j-1,k)
                               - state(i  ,j-1,k,n)*mac(i  ,j-1,k))
           // divu = 0
           //- dt/(3.e0)   *s(i,j-1,k,n)*divu_cc(i,j-1,k)
            + dt/(3.e0*dx)*s(i,j-1,k,n)*(mac(i+1,j-1,k)-mac(i,j-1,k));

        hi1  = hi - dt/(3.e0*dx)*(state(i+1,j,k,n)*mac(i+1,j,k)
                                - state(i  ,j,k,n)*mac(i  ,j,k))
           // divu = 0
           //- dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
            + dt/(3.e0*dx)*s(i,j,k,n)*(mac(i+1,j,k) - mac(i,j,k));
    } else {
        lo1 = lo - dt/(6.e0*dx)*(mac  (i+1,j-1,k)   +   mac(i,j-1,k))
                              * (state(i+1,j-1,k,n) - state(i,j-1,k,n));
        hi1 = hi - dt/(6.e0*dx)*(mac  (i+1,j  ,k)   +   mac(i,j  ,k))
                              * (state(i+1,j  ,k,n) - state(i,j  ,k,n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_corner_couple_zx (amrex::Real& lo1, amrex::Real& hi1,
                               int i, int j, int k, int n, amrex::Real dt, amrex::Real dx,
                               bool iconserv, amrex::Real lo, amrex::Real hi,
                               amrex::Array4<amrex::Real const> const& s,
                               amrex::Array4<amrex::Real const> const& mac,
                               amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo - dt/(3.e0*dx)*(state(i+1,j,k-1,n)*mac(i+1,j,k-1)
                               - state(i  ,j,k-1,n)*mac(i  ,j,k-1))
            // divu = 0
            //- dt/(3.e0)   *s(i,j,k-1,n)*divu_cc(i,j,k-1)
            + dt/(3.e0*dx)*s(i,j,k-1,n)*(mac(i+1,j,k-1)-mac(i,j,k-1));

        hi1 = hi - dt/(3.e0*dx)*(state(i+1,j,k,n)*mac(i+1,j,k)
                               - state(i  ,j,k,n)*mac(i  ,j,k))
            // divu = 0
            //- dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
            + dt/(3.e0*dx)*s(i,j,k,n)*(mac(i+1,j,k) - mac(i,j,k));
    } else {
        lo1 = lo - dt/(6.e0*dx)*(  mac(i+1,j,k-1)   +   mac(i,j,k-1))
                              * (state(i+1,j,k-1,n) - state(i,j,k-1,n));
        hi1 = hi - dt/(6.e0*dx)*(  mac(i+1,j,k  )   +   mac(i,j,k))
                              * (state(i+1,j,k  ,n) - state(i,j,k  ,n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_corner_couple_xy (amrex::Real& lo1, amrex::Real& hi1,
                               int i, int j, int k, int n, amrex::Real dt, amrex::Real dy,
                               bool iconserv, amrex::Real lo, amrex::Real hi,
                               amrex::Array4<amrex::Real const> const& s,
                               amrex::Array4<amrex::Real const> const& mac,
                               amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo - dt/(3.e0*dy)*(state(i-1,j+1,k,n)*mac(i-1,j+1,k)
                               - state(i-1,j  ,k,n)*mac(i-1,j  ,k))
            // divu = 0
            //- dt/3.e0     *s(i-1,j,k,n)*divu_cc(i-1,j,k)
            + dt/(3.e0*dy)*s(i-1,j,k,n)*(mac(i-1,j+1,k)-mac(i-1,j,k));

        hi1 = hi - dt/(3.e0*dy)*(state(i,j+1,k,n)*mac(i,j+1,k)
                               - state(i,j  ,k,n)*mac(i,j  ,k))
            // divu = 0
            //- dt/3.e0     *s(i,j,k,n)*divu_cc(i,j,k)
            + dt/(3.e0*dy)*s(i,j,k,n)*(mac(i,j+1,k)-mac(i,j,k));
    } else {
        lo1 = lo - dt/(6.e0*dy)*(  mac(i-1,j+1,k)   +   mac(i-1,j,k))
                              * (state(i-1,j+1,k,n) - state(i-1,j,k,n));
        hi1 = hi - dt/(6.e0*dy)*(mac  (i  ,j+1,k)   +   mac(i  ,j,k))
                              * (state(i  ,j+1,k,n) - state(i  ,j,k,n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_corner_couple_zy (amrex::Real& lo1, amrex::Real& hi1,
                               int i, int j, int k, int n, amrex::Real dt, amrex::Real dy,
                               bool iconserv, amrex::Real lo, amrex::Real hi,
                               amrex::Array4<amrex::Real const> const& s,
                               amrex::Array4<amrex::Real const> const& mac,
                               amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo - dt/(3.e0*dy)*(state(i,j+1,k-1,n)*mac(i,j+1,k-1)
                               - state(i,j  ,k-1,n)*mac(i,j  ,k-1))
            // divu = 0
            //- dt/(3.e0)   *s(i,j,k-1,n)*divu_cc(i,j,k-1)
            + dt/(3.e0*dy)*s(i,j,k-1,n)*(mac(i,j+1,k-1)-mac(i,j,k-1));

        hi1 = hi - dt/(3.e0*dy)*(state(i,j+1,k,n)*mac(i,j+1,k)
                               - state(i,j  ,k,n)*mac(i,j  ,k))
            // divu = 0
            //- dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
            + dt/(3.e0*dy)*s(i,j,k,n)*(mac(i,j+1,k) - mac(i,j,k));
    } else {
        lo1 = lo - dt/(6.e0*dy)*(  mac(i,j+1,k-1)   +   mac(i,j,k-1))
                              * (state(i,j+1,k-1,n) - state(i,j,k-1,n));
        hi1 = hi - dt/(6.e0*dy)*(  mac(i,j+1,k)     +   mac(i,j,k))
                              * (state(i,j+1,k  ,n) - state(i,j,k  ,n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_corner_couple_xz (amrex::Real& lo1, amrex::Real& hi1,
                               int i, int j, int k, int n, amrex::Real dt, amrex::Real dz,
                               bool iconserv, amrex::Real lo, amrex::Real hi,
                               amrex::Array4<amrex::Real const> const& s,
                               amrex::Array4<amrex::Real const> const& mac,
                               amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo - dt/(3.e0*dz)*(state(i-1,j,k+1,n)*mac(i-1,j,k+1)
                               - state(i-1,j,k  ,n)*mac(i-1,j,k))
            // divu = 0
            //- dt/(3.e0)   *s(i-1,j,k,n)*divu_cc(i-1,j,k)
            + dt/(3.e0*dz)*s(i-1,j,k,n)*(mac(i-1,j,k+1) - mac(i-1,j,k));

        hi1 = hi - dt/(3.e0*dz)*(state(i,j,k+1,n)*mac(i,j,k+1)
                               - state(i,j,k  ,n)*mac(i,j,k))
            // divu = 0
            //- dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
            + dt/(3.e0*dz)*s(i,j,k,n)*(mac(i,j,k+1) - mac(i,j,k));
    } else {
        lo1 = lo - dt/(6.e0*dz)*(  mac(i-1,j,k+1)   +   mac(i-1,j,k))
                              * (state(i-1,j,k+1,n) - state(i-1,j,k,n));
        hi1 = hi - dt/(6.e0*dz)*(  mac(i  ,j,k+1)   +   mac(i  ,j,k))
                              * (state(i  ,j,k+1,n) - state(i  ,j,k,n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_corner_couple_yz (amrex::Real& lo1, amrex::Real& hi1,
                               int i, int j, int k, int n, amrex::Real dt, amrex::Real dz,
                               bool iconserv, amrex::Real lo, amrex::Real hi,
                               amrex::Array4<amrex::Real const> const& s,
                               amrex::Array4<amrex::Real const> const& mac,
                               amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo - dt/(3.e0*dz)*(state(i,j-1,k+1,n)*mac(i,j-1,k+1)
                               - state(i,j-1,k  ,n)*mac(i,j-1,k))
            // divu = 0
            //- dt/(3.e0)   *s(i,j-1,k,n)*divu_cc(i,j-1,k)
            + dt/(3.e0*dz)*s(i,j-1,k,n)*(mac(i,j-1,k+1)-mac(i,j-1,k));

        hi1 = hi - dt/(3.e0*dz)*(state(i,j,k+1,n)*mac(i,j,k+1)
                               - state(i,j,k  ,n)*mac(i,j,k))
            // divu = 0
            //- dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
            + dt/(3.e0*dz)*s(i,j,k,n)*(mac(i,j,k+1) - mac(i,j,k));
    } else {
        lo1 = lo - dt/(6.e0*dz)*(  mac(i,j-1,k+1)   +   mac(i,j-1,k))
                              * (state(i,j-1,k+1,n) - state(i,j-1,k,n));
        hi1 = hi - dt/(6.e0*dz)*(  mac(i,j  ,k+1)   +   mac(i,j  ,k))
                              * (state(i,j  ,k+1,n) - state(i,j  ,k,n));
    }
}

}
#endif 
