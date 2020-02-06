#ifndef GODUNOV_ADVECTION_H
#define GODUNOV_ADVECTION_H
#include <incflo.H> 
#include <AMReX_Gpu.H> 
#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H> 
#include <AMReX_BC_TYPES.H>
#include <AMReX_Array.H>
#include <iomanip>

/* This header file contains the inlined __host__ __device__ functions required for 
   the scalar advection routines for 3D Godunov. It also contains fuction declarations 
   for controlling host functions.  */


namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real vanLeer(const amrex::Real a, const amrex::Real b, const amrex::Real c)
{   
    amrex::Real dsc = 0.5e0*(b - c);
    amrex::Real dsl = 2.e0*(a - c);
    amrex::Real dsr = 2.e0*(b - a);
    return (dsl*dsr > 0.)? std::copysign(1., dsc)*amrex::min(std::abs(dsc),
                                 amrex::min(std::abs(dsl), std::abs(dsr))) : 0.;
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

        if(bclo == BCType::ext_dir || bclo == BCType::hoextrap){
            if(i == domlo){
                 sm = s(domlo-1,j,k,n);
                 sedge2 = -0.2*s(domlo-1, j, k, n) + 0.75*s(domlo,j, k, n) 
                       +  0.5*s(domlo+1,j,k,n) - 0.05*s(domlo+2,j,k,n);
                sedge2 = amrex::max(sedge2, amrex::min(s(domlo+1,j,k,n), s(domlo,j,k,n)));
                sedge2 = amrex::min(sedge2, amrex::max(s(domlo+1,j,k,n), s(domlo,j,k,n)));
                sp = sedge2;
            }    
            else if(i == domlo+1){
                sedge1 = -0.2*s(domlo-1, j, k, n) + 0.75*s(domlo,j,k, n) 
                       +  0.5*s(domlo+1,j,k,n) - 0.05*s(domlo+2,j,k,n);
                sedge1 = amrex::max(sedge1, amrex::min(s(domlo+1,j,k,n), s(domlo,j,k,n)));
                sedge1 = amrex::min(sedge1, amrex::max(s(domlo+1,j,k,n), s(domlo,j,k,n)));

                sp = sedge2;
                sm = sedge1;
                if( (sp - s(domlo+1,j,k,n))*(s(domlo+1,j,k,n) - sm) <= 0.){
                    sp = s(domlo+1,j,k,n);
                    sm = s(domlo+1,j,k,n);
                }
                else if(std::abs(sp - s(domlo+1,j,k,n)) >= 2.*std::abs(sm - s(domlo+1,j,k,n))){
                    sp = 3.*s(domlo+1,j,k,n) - 2.*sm;
                }
                else if(std::abs(sm - s(domlo+1,j,k,n)) >= 2.*std::abs(sp - s(domlo+1,j,k,n))){
                    sm = 3.*s(domlo+1,j,k,n) - 2.*sp;
                } 
            }

        }

      if(bchi == BCType::ext_dir || bchi == BCType::hoextrap){
          if(i == domhi){ 
            sp = s(domhi,j,k+1, n);
            sedge1 =  -0.2*s(domhi+1, j, k, n) + 0.75*s(domhi,j,k, n) 
                       +  0.5*s(domhi-1,j,k,n) - 0.05*s(domhi-2,j,k,n);
            sedge1 = amrex::max(sedge1, amrex::min(s(domhi-1,j,k,n), s(domhi,j,k,n)));
            sedge1 = amrex::min(sedge1, amrex::max(s(domhi-1,j,k,n), s(domhi,j,k,n)));
            sm = sedge1;
          }
          else if(i == domhi - 1){
                sedge2 =  -0.2*s(domhi+1, j,k, n) + 0.75*s(domhi,j, k, n) 
                       +  0.5*s(domhi-1,j,k,n) - 0.05*s(domhi-2,j,k,n);
                sedge2 = amrex::max(sedge2, amrex::min(s(domhi-1,j,k,n), s(domhi,j,k,n)));
                sedge2 = amrex::min(sedge2, amrex::max(s(domhi-1,j,k,n), s(domhi,j,k,n)));
                sp = sedge2;
                sm = sedge1;
                if( (sp - s(domhi-1,j,k,n))*(s(domhi-1,j,k,n) - sm) <= 0.){
                    sp = s(domhi-1,j,k,n);
                    sm = s(domhi-1,j,k,n);
                }
                else if(std::abs(sp - s(domhi-1,j,k,n)) >= 2.*std::abs(sm - s(domhi-1,j,k,n))){
                    sp = 3.*s(domhi-1,j,k,n) - 2.*sm;
                }
                else if(std::abs(sm - s(domhi-1,j,k,n)) >= 2.*std::abs(sp - s(domhi-1,j,k,n))){
                    sm = 3.*s(domhi-1,j,k,n) - 2.*sp;
                } 
    
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

        if(bclo == BCType::ext_dir || bclo == BCType::hoextrap){
            if(j == domlo){
                 sm = s(i,domlo-1,k,n);
                 sedge2 = -0.2*s(i, domlo-1, k, n) + 0.75*s(i,domlo, k, n) 
                       +  0.5*s(i,domlo+1,k,n) - 0.05*s(i,domlo+2,k,n);
                sedge2 = amrex::max(sedge2, amrex::min(s(i,domlo+1,k,n), s(i,domlo,k,n)));
                sedge2 = amrex::min(sedge2, amrex::max(s(i,domlo+1,k,n), s(i,domlo,k,n)));
                sp = sedge2;
            }    
            else if(j == domlo+1){
                sedge1 = -0.2*s(i, domlo-1, k, n) + 0.75*s(i,domlo, k, n) 
                       +  0.5*s(i,domlo+1,k,n) - 0.05*s(i,domlo+2,k,n);
                sedge1 = amrex::max(sedge1, amrex::min(s(i,domlo+1,k,n), s(i,domlo,k,n)));
                sedge1 = amrex::min(sedge1, amrex::max(s(i,domlo+1,k,n), s(i,domlo,k,n)));

                sp = sedge2;
                sm = sedge1;
                if( (sp - s(i,domlo+1,k,n))*(s(i,domlo+1,k,n) - sm) <= 0.){
                    sp = s(i,domlo+1,k,n);
                    sm = s(i,domlo+1,k,n);
                }
                else if(std::abs(sp - s(i,domlo+1,k,n)) >= 2.*std::abs(sm - s(i,domlo+1,k,n))){
                    sp = 3.*s(i,domlo+1,k,n) - 2.*sm;
                }
                else if(std::abs(sm - s(i,domlo+1,k,n)) >= 2.*std::abs(sp - s(i,domlo+1,k,n))){
                    sm = 3.*s(i,domlo+1,k,n) - 2.*sp;
                } 
            }

        }

      if(bchi == BCType::ext_dir || bchi == BCType::hoextrap){
          if(j == domhi){ 
            sp = s(i,domhi+1,k, n);
            sedge1 =  -0.2*s(i, domhi+1,k, n) + 0.75*s(i,domhi, k, n) 
                       +  0.5*s(i,domhi-1,k,n) - 0.05*s(i,domhi-2,k,n);
            sedge1 = amrex::max(sedge1, amrex::min(s(i,domhi-1,k,n), s(i,domhi,k,n)));
            sedge1 = amrex::min(sedge1, amrex::max(s(i,domhi-1,k,n), s(i,domhi,k,n)));
            sm = sedge1;
          }
          else if(j == domhi - 1){
            sedge2 =  -0.2*s(i, domhi+1,k, n) + 0.75*s(i,domhi, k, n) 
                       +  0.5*s(i,domhi-1,k,n) - 0.05*s(i,domhi-2,k,n);
            sedge2 = amrex::max(sedge2, amrex::min(s(i,domhi-1,k,n), s(i,domhi,k,n)));
            sedge2 = amrex::min(sedge2, amrex::max(s(i,domhi-1,k,n), s(i,domhi,k,n)));
                sp = sedge2;
                sm = sedge1;
                if( (sp - s(i,domhi-1,k,n))*(s(i,domhi-1,k,n) - sm) <= 0.){
                    sp = s(i,domhi-1,k,n);
                    sm = s(i,domhi-1,k,n);
                }
                else if(std::abs(sp - s(i,domhi-1,k,n)) >= 2.*std::abs(sm - s(i,domhi-1,k,n))){
                    sp = 3.*s(i,domhi-1,k,n) - 2.*sm;
                }
                else if(std::abs(sm - s(i,domhi-1,k,n)) >= 2.*std::abs(sp - s(i,domhi-1,k,n))){
                    sm = 3.*s(i,domhi-1,k,n) - 2.*sp;
                } 
    
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

        if(bclo == BCType::ext_dir || bclo == BCType::hoextrap){
            if(k == domlo){
                 sm = s(i,j,domlo-1,n);
                 sedge2 = -0.2*s(i, j, domlo-1, n) + 0.75*s(i,j, domlo, n) 
                       +  0.5*s(i,j,domlo+1,n) - 0.05*s(i,j,domlo+2,n);
                sedge2 = amrex::max(sedge2, amrex::min(s(i,j,domlo+1,n), s(i,j,domlo,n)));
                sedge2 = amrex::min(sedge2, amrex::max(s(i,j,domlo+1,n), s(i,j,domlo,n)));
                sp = sedge2;
            }    
            else if(k == domlo+1){
                sedge1 = -0.2*s(i, j, domlo-1, n) + 0.75*s(i,j, domlo, n) 
                       +  0.5*s(i,j,domlo+1,n) - 0.05*s(i,j,domlo+2,n);
                sedge1 = amrex::max(sedge1, amrex::min(s(i,j,domlo+1,n), s(i,j,domlo,n)));
                sedge1 = amrex::min(sedge1, amrex::max(s(i,j,domlo+1,n), s(i,j,domlo,n)));

                sp = sedge2;
                sm = sedge1;
                if( (sp - s(i,j,domlo+1,n))*(s(i,j,domlo+1,n) - sm) <= 0.){
                    sp = s(i,j,domlo+1,n);
                    sm = s(i,j,domlo+1,n);
                }
                else if(std::abs(sp - s(i,j,domlo+1,n)) >= 2.*std::abs(sm - s(i,j,domlo+1,n))){
                    sp = 3.*s(i,j,domlo+1,n) - 2.*sm;
                }
                else if(std::abs(sm - s(i,j,domlo+1,n)) >= 2.*std::abs(sp - s(i,j,domlo+1,n))){
                    sm = 3.*s(i,j,domlo+1,n) - 2.*sp;
                } 
            }

        }

      if(bchi == BCType::ext_dir || bchi == BCType::hoextrap){
          if(k == domhi){ 
            sp = s(i,j,domhi+1, n);
            sedge1 =  -0.2*s(i, j, domhi+1, n) + 0.75*s(i,j, domhi, n) 
                       +  0.5*s(i,j,domhi-1,n) - 0.05*s(i,j,domhi-2,n);
            sedge1 = amrex::max(sedge1, amrex::min(s(i,j,domhi-1,n), s(i,j,domhi,n)));
            sedge1 = amrex::min(sedge1, amrex::max(s(i,j,domhi-1,n), s(i,j,domhi,n)));
            sm = sedge1;
          }
          else if(k == domhi - 1){
                sedge2 =  -0.2*s(i, j, domhi+1, n) + 0.75*s(i,j, domhi, n) 
                       +  0.5*s(i,j,domhi-1,n) - 0.05*s(i,j,domhi-2,n);
                sedge2 = amrex::max(sedge2, amrex::min(s(i,j,domhi-1,n), s(i,j,domhi,n)));
                sedge2 = amrex::min(sedge2, amrex::max(s(i,j,domhi-1,n), s(i,j,domhi,n)));
                sp = sedge2;
                sm = sedge1;
                if( (sp - s(i,j,domhi-1,n))*(s(i,j,domhi-1,n) - sm) <= 0.){
                    sp = s(i,j,domhi-1,n);
                    sm = s(i,j,domhi-1,n);
                }
                else if(std::abs(sp - s(i,j,domhi-1,n)) >= 2.*std::abs(sm - s(i,j,domhi-1,n))){
                    sp = 3.*s(i,j,domhi-1,n) - 2.*sm;
                }
                else if(std::abs(sm - s(i,j,domhi-1,n)) >= 2.*std::abs(sp - s(i,j,domhi-1,n))){
                    sm = 3.*s(i,j,domhi-1,n) - 2.*sp;
                } 
    
           }
      }
} 

// Right now only ppm type 1 is supported on GPU
// This version is called before the MAC projection, when we use the cell-centered velocity
//      for upwinding
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_ppm_pred_x(const int i, const int j, const int k, const int n,
                             const amrex::Real dt, const amrex::Real dx, 
                             const amrex::Array4<const amrex::Real> &S, 
                             const amrex::Array4<const amrex::Real> &vel, 
                             const amrex::Array4<amrex::Real> &Im, 
                             const amrex::Array4<amrex::Real> &Ip,
                             const amrex::BCRec bc, const int domlo, const int domhi)
{
    using namespace amrex;

    amrex::Real sedge1, sedge2, sp, sm;
    amrex::Real s1 = S(i,j,k,n);
    constexpr amrex::Real eps = 1.e-6;

    {
        amrex::Real s2 = S(i-1,j,k,n);
        amrex::Real d1 = vanLeer(S(i,j,k,n), S(i+1,j,k,n), S(i-1,j,k,n));
        amrex::Real d2 = vanLeer(S(i-1,j,k,n), S(i,j,k,n), S(i-2,j,k,n));
        sedge1 = 0.5e0*(S(i,j,k,n) + S(i-1,j,k,n)) - 1.e0/6.e0*(d1 - d2);
        sedge1 = amrex::max(sedge1, amrex::min(S(i,j,k,n), S(i-1,j,k,n)));
        sedge1 = amrex::min(sedge1, amrex::max(S(i,j,k,n), S(i-1,j,k,n)));

        s2 = s1;
        s1 = S(i+1,j,k,n);

        d2 = d1;
        d1 = vanLeer(S(i+1,j,k,n), S(i+2, j, k,n), S(i,j,k,n));
        sedge2 = 0.5e0*(S(i+1,j,k,n) + S(i,j,k,n)) - 1.e0/6.e0*(d1 - d2);
        sedge2 = amrex::max(sedge2, amrex::min(S(i+1,j,k,n), S(i,j,k,n)));
        sedge2 = amrex::min(sedge2, amrex::max(S(i+1,j,k,n), S(i,j,k,n)));
        
        sm = sedge1;
        sp = sedge2;
        if ((sedge2-S(i,j,k,n))*(S(i,j,k,n)-sedge1) <= 0.e0){
            sp = S(i,j,k,n);
            sm = S(i,j,k,n);
        }
        else if (std::abs(sp-S(i,j,k,n)) >= 2.e0*std::abs(sm-S(i,j,k,n)))
            sp = 3.e0*S(i,j,k,n) - 2.e0*sedge1;
        else if (std::abs(sm-S(i,j,k,n)) >=  2.e0*std::abs(sp-S(i,j,k,n))) 
            sm = 3.e0*S(i,j,k,n) - 2.e0*sedge2;
        
        Godunov_ppm_xbc(i, j, k, n, sm, sp, sedge1,sedge2, S, bc.lo(0), bc.hi(0), domlo, domhi);
        amrex::Real s6 = 6.e0*s2- 3.e0*(sm + sp);
        amrex::Real sigma = std::abs(vel(i,j,k,0))*dt/dx;

        if (vel(i,j,k,0) > eps)
        {
            Ip(i,j,k,n) = sp - (0.5*sigma)*((sp - sm) 
                             - (1.e0 -2.e0/3.e0*sigma)*s6);
            Im(i,j,k,n) = S(i,j,k,n);
        }
        else
        {
            Ip(i,j,k,n) = S(i,j,k,n);
            Im(i,j,k,n) = sm + (0.5*sigma)*((sp-sm)
                             + (1.e0 - 2.e0/3.e0*sigma)*s6);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_ppm_pred_y(const int i, const int j, const int k, const int n,
                             const amrex::Real dt, const amrex::Real dx,
                             const amrex::Array4<const amrex::Real> &S,
                             const amrex::Array4<const amrex::Real> &vel,
                             const amrex::Array4<amrex::Real> &Im,
                             const amrex::Array4<amrex::Real> &Ip,
                             const amrex::BCRec bc, const int domlo, const int domhi)
{
    using namespace amrex;

    amrex::Real sedge1, sedge2, sp, sm;
    amrex::Real s1 = S(i,j,k,n);
    constexpr amrex::Real eps = 1.e-6;

    {
        amrex::Real s2 = S(i,j-1,k,n);
        amrex::Real d1 = vanLeer(s1, S(i,j+1,k,n), s2);
        amrex::Real d2 = vanLeer(s2, s1, S(i,j-2,k,n));
        sedge1 = 0.5e0*(s1 + s2) - 1.e0/6.e0*(d1 - d2);
        sedge1 = amrex::max(sedge1, amrex::min(s1, s2));
        sedge1 = amrex::min(sedge1, amrex::max(s1, s2));

        amrex::Real s3 = s2;
        s2 = s1;
        s1 = S(i,j+1,k,n);
        
        d1 = vanLeer(s1, S(i, j+2, k,n), S(i,j,k,n));
        d2 = vanLeer(S(i,j,k,n), s1, s3);
        sedge2 = 0.5e0*(s1 + S(i,j,k,n)) - 1.e0/6.e0*(d1 - d2);
        sedge2 = amrex::max(sedge2, amrex::min(s1, S(i,j,k,n)));
        sedge2 = amrex::min(sedge2, amrex::max(s1, S(i,j,k,n)));
        
        sm = sedge1;
        sp = sedge2;
        if ((sedge2-S(i,j,k,n))*(S(i,j,k,n)-sedge1) < 0.e0){
            sp = S(i,j,k,n);
            sm = S(i,j,k,n);
        }
        else if (std::abs(sedge2-S(i,j,k,n)) >= 2.e0*std::abs(sedge1-S(i,j,k,n)))
            sp = 3.e0*S(i,j,k,n) - 2.e0*sedge1;
        else if (std::abs(sedge1-S(i,j,k,n)) >=  2.e0*std::abs(sedge2-S(i,j,k,n))) 
            sm = 3.e0*S(i,j,k,n) - 2.e0*sedge2;
        Godunov_ppm_ybc(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(1), bc.hi(1), domlo, domhi);

        amrex::Real s6 = 6.e0*S(i,j,k,n)- 3.e0*(sm + sp);

        amrex::Real sigma = std::abs(vel(i,j,k,1))*dt/dx;

        if(vel(i,j,k,1) > eps)
        {
            Ip(i,j,k,n) = sp - (0.5*sigma)*((sp - sm) 
                             - (1.e0 -2.e0/3.e0*sigma)*s6);
            Im(i,j,k,n) = S(i,j,k,n);
        }
        else
        {
            Im(i,j,k,n) = sm + (0.5*sigma)*((sp-sm) 
                             + (1.e0 - 2.e0/3.e0*sigma)*s6);
            Ip(i,j,k,n) = S(i,j,k,n);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_ppm_pred_z(const int i, const int j, const int k, const int n,
                             const amrex::Real dt, const amrex::Real dx,
                             const amrex::Array4<const amrex::Real> &S,
                             const amrex::Array4<const amrex::Real> &vel,
                             const amrex::Array4<amrex::Real> &Im,
                             const amrex::Array4<amrex::Real> &Ip,
                             const amrex::BCRec bc, const int domlo, const int domhi)
{
    using namespace amrex;

    amrex::Real sedge1, sedge2, sp, sm;
    amrex::Real s1 = S(i,j,k,n);
    constexpr amrex::Real eps = 1.e-6;

    {
        amrex::Real s2 = S(i,j,k-1,n);

        amrex::Real d1 = vanLeer(S(i,j,k,n), S(i,j,k+1,n), S(i,j,k-1,n));
        amrex::Real d2 = vanLeer(S(i,j,k-1,n), S(i,j,k,n), S(i,j,k-2,n));

        sedge1 = 0.5e0*(S(i,j,k,n) + S(i,j,k-1,n)) - 1.e0/6.e0*(d1 - d2);
        sedge1 = amrex::max(sedge1, amrex::min(S(i,j,k,n), S(i,j,k-1,n)));
        sedge1 = amrex::min(sedge1, amrex::max(S(i,j,k,n), S(i,j,k-1,n)));

        s2 = s1;
        s1 = S(i,j,k+1,n);

        d2 = d1;
        d1 = vanLeer(S(i,j,k+1,n), S(i, j, k+2,n), S(i,j,k,n));

        sedge2 = 0.5e0*(S(i,j,k+1,n) + S(i,j,k,n)) - 1.e0/6.e0*(d1 - d2);
        sedge2 = amrex::max(sedge2, amrex::min(S(i,j,k+1,n), S(i,j,k,n)));
        sedge2 = amrex::min(sedge2, amrex::max(S(i,j,k+1,n), S(i,j,k,n)));
        
        sm = sedge1;
        sp = sedge2;
        if ((sp-S(i,j,k,n))*(S(i,j,k,n)-sm) <= 0.e0){
            sp = S(i,j,k,n);
            sm = S(i,j,k,n);
        }
        else if (std::abs(sp-S(i,j,k,n)) >= 2.e0*std::abs(sm-S(i,j,k,n)))
            sp = 3.e0*S(i,j,k,n) - 2.e0*sm;
        else if (std::abs(sm-S(i,j,k,n)) >=  2.e0*std::abs(sp-S(i,j,k,n))) 
            sm = 3.e0*S(i,j,k,n) - 2.e0*sp;

        Godunov_ppm_zbc(i, j, k, n, sm, sp, sedge1, sedge2,S, bc.lo(2), bc.hi(2), domlo, domhi);

        amrex::Real s6 = 6.e0*S(i,j,k,n)- 3.e0*(sm + sp);
        amrex::Real sigma = std::abs(vel(i,j,k,2))*dt/dx;

        if(vel(i,j,k,2) > eps)
        { 
            Ip(i,j,k,n) = sp - (0.5*sigma)*((sp - sm)
                             - (1.e0 -2.e0/3.e0*sigma)*s6);
            Im(i,j,k,n) = s2;
        } 
        else
        { 
            Im(i,j,k,n) = sm + (0.5*sigma)*((sp-sm) 
                             + (1.e0 - 2.e0/3.e0*sigma)*s6);
            Ip(i,j,k,n) = s2;
        } 
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
    using namespace amrex;

    amrex::Real sedge1, sedge2, sp, sm;
    amrex::Real s1 = S(i,j,k,n);
    
    {
        amrex::Real s2 = S(i-1,j,k,n);
        amrex::Real d1 = vanLeer(S(i,j,k,n), S(i+1,j,k,n), S(i-1,j,k,n));
        amrex::Real d2 = vanLeer(S(i-1,j,k,n), S(i,j,k,n), S(i-2,j,k,n));
        sedge1 = 0.5e0*(S(i,j,k,n) + S(i-1,j,k,n)) - 1.e0/6.e0*(d1 - d2);
        sedge1 = amrex::max(sedge1, amrex::min(S(i,j,k,n), S(i-1,j,k,n)));
        sedge1 = amrex::min(sedge1, amrex::max(S(i,j,k,n), S(i-1,j,k,n)));

        s2 = s1;
        s1 = S(i+1,j,k,n);

        d2 = d1;
        d1 = vanLeer(S(i+1,j,k,n), S(i+2, j, k,n), S(i,j,k,n));
        sedge2 = 0.5e0*(S(i+1,j,k,n) + S(i,j,k,n)) - 1.e0/6.e0*(d1 - d2);
        sedge2 = amrex::max(sedge2, amrex::min(S(i+1,j,k,n), S(i,j,k,n)));
        sedge2 = amrex::min(sedge2, amrex::max(S(i+1,j,k,n), S(i,j,k,n)));
        
        sm = sedge1;
        sp = sedge2;
        if ((sedge2-S(i,j,k,n))*(S(i,j,k,n)-sedge1) <= 0.e0){
            sp = S(i,j,k,n);
            sm = S(i,j,k,n);
        }
        else if (std::abs(sp-S(i,j,k,n)) >= 2.e0*std::abs(sm-S(i,j,k,n)))
            sp = 3.e0*S(i,j,k,n) - 2.e0*sedge1;
        else if (std::abs(sm-S(i,j,k,n)) >=  2.e0*std::abs(sp-S(i,j,k,n))) 
            sm = 3.e0*S(i,j,k,n) - 2.e0*sedge2;
        
        Godunov_ppm_xbc(i, j, k, n, sm, sp, sedge1,sedge2, S, bc.lo(0), bc.hi(0), domlo, domhi);
        amrex::Real s6 = 6.e0*s2- 3.e0*(sm + sp);
        amrex::Real sigmap = std::abs(vel_edge(i+1,j,k))*dt/dx;
        amrex::Real sigmam = std::abs(vel_edge(i,j,k))*dt/dx;

        if(vel_edge(i+1,j,k) > 1.e-06){
            Ip = sp - (0.5*sigmap)*((sp - sm)
                             - (1.e0 -2.e0/3.e0*sigmap)*s6);
        }
        else
            Ip = S(i,j,k,n);

        if(vel_edge(i,j,k) < -1.e-06)
            Im = sm + (0.5*sigmam)*((sp-sm)
                             + (1.e0 - 2.e0/3.e0*sigmam)*s6);
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
    using namespace amrex;

    amrex::Real sedge1, sedge2, sp, sm;
    amrex::Real s1 = S(i,j,k,n);

    {
        amrex::Real s2 = S(i,j-1,k,n);
        amrex::Real d1 = vanLeer(s1, S(i,j+1,k,n), s2);
        amrex::Real d2 = vanLeer(s2, s1, S(i,j-2,k,n));
        sedge1 = 0.5e0*(s1 + s2) - 1.e0/6.e0*(d1 - d2);
        sedge1 = amrex::max(sedge1, amrex::min(s1, s2));
        sedge1 = amrex::min(sedge1, amrex::max(s1, s2));

        amrex::Real s3 = s2;
        s2 = s1;
        s1 = S(i,j+1,k,n);
        
        d1 = vanLeer(s1, S(i, j+2, k,n), S(i,j,k,n));
        d2 = vanLeer(S(i,j,k,n), s1, s3);
        sedge2 = 0.5e0*(s1 + S(i,j,k,n)) - 1.e0/6.e0*(d1 - d2);
        sedge2 = amrex::max(sedge2, amrex::min(s1, S(i,j,k,n)));
        sedge2 = amrex::min(sedge2, amrex::max(s1, S(i,j,k,n)));
        
        sm = sedge1;
        sp = sedge2;
        if ((sedge2-S(i,j,k,n))*(S(i,j,k,n)-sedge1) < 0.e0){
            sp = S(i,j,k,n);
            sm = S(i,j,k,n);
        }
        else if (std::abs(sedge2-S(i,j,k,n)) >= 2.e0*std::abs(sedge1-S(i,j,k,n)))
            sp = 3.e0*S(i,j,k,n) - 2.e0*sedge1;
        else if (std::abs(sedge1-S(i,j,k,n)) >=  2.e0*std::abs(sedge2-S(i,j,k,n))) 
            sm = 3.e0*S(i,j,k,n) - 2.e0*sedge2;
        Godunov_ppm_ybc(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(1), bc.hi(1), domlo, domhi);

        amrex::Real s6 = 6.e0*S(i,j,k,n)- 3.e0*(sm + sp);
        amrex::Real sigmap = std::abs(vel_edge(i,j+1,k))*dt/dx;
        amrex::Real sigmam = std::abs(vel_edge(i,j,k))*dt/dx;
        if(vel_edge(i,j+1,k) > 1.e-06)
            Ip = sp - (0.5*sigmap)*((sp - sm)
                             - (1.e0 -2.e0/3.e0*sigmap)*s6);
        else
            Ip = S(i,j,k,n);

        if(vel_edge(i,j,k) < -1.e-06)
            Im = sm + (0.5*sigmam)*((sp-sm)
                             + (1.e0 - 2.e0/3.e0*sigmam)*s6);
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
    using namespace amrex;

    amrex::Real sedge1, sedge2, sp, sm;
    amrex::Real s1 = S(i,j,k,n);

    {
        amrex::Real s2 = S(i,j,k-1,n);

        amrex::Real d1 = vanLeer(S(i,j,k,n), S(i,j,k+1,n), S(i,j,k-1,n));
        amrex::Real d2 = vanLeer(S(i,j,k-1,n), S(i,j,k,n), S(i,j,k-2,n));

        sedge1 = 0.5e0*(S(i,j,k,n) + S(i,j,k-1,n)) - 1.e0/6.e0*(d1 - d2);
        sedge1 = amrex::max(sedge1, amrex::min(S(i,j,k,n), S(i,j,k-1,n)));
        sedge1 = amrex::min(sedge1, amrex::max(S(i,j,k,n), S(i,j,k-1,n)));

        s2 = s1;
        s1 = S(i,j,k+1,n);

        d2 = d1;
        d1 = vanLeer(S(i,j,k+1,n), S(i, j, k+2,n), S(i,j,k,n));

        sedge2 = 0.5e0*(S(i,j,k+1,n) + S(i,j,k,n)) - 1.e0/6.e0*(d1 - d2);
        sedge2 = amrex::max(sedge2, amrex::min(S(i,j,k+1,n), S(i,j,k,n)));
        sedge2 = amrex::min(sedge2, amrex::max(S(i,j,k+1,n), S(i,j,k,n)));
        
        sm = sedge1;
        sp = sedge2;
        if ((sp-S(i,j,k,n))*(S(i,j,k,n)-sm) <= 0.e0){
            sp = S(i,j,k,n);
            sm = S(i,j,k,n);
        }
        else if (std::abs(sp-S(i,j,k,n)) >= 2.e0*std::abs(sm-S(i,j,k,n)))
            sp = 3.e0*S(i,j,k,n) - 2.e0*sm;
        else if (std::abs(sm-S(i,j,k,n)) >=  2.e0*std::abs(sp-S(i,j,k,n))) 
            sm = 3.e0*S(i,j,k,n) - 2.e0*sp;

        Godunov_ppm_zbc(i, j, k, n, sm, sp, sedge1, sedge2,S, bc.lo(2), bc.hi(2), domlo, domhi);

        amrex::Real s6 = 6.e0*S(i,j,k,n)- 3.e0*(sm + sp);
        amrex::Real sigmap = std::abs(vel_edge(i,j,k+1))*dt/dx;
        amrex::Real sigmam = std::abs(vel_edge(i,j,k))*dt/dx;
        if(vel_edge(i,j,k+1) > 1.e-06)
            Ip = sp - (0.5*sigmap)*((sp - sm)
                             - (1.e0 -2.e0/3.e0*sigmap)*s6);
        else
            Ip = s2;


        if(vel_edge(i,j,k) < -1.e-06)
            Im = sm + (0.5*sigmam)*((sp-sm)
                             + (1.e0 - 2.e0/3.e0*sigmam)*s6);
        else
            Im = s2;
    }
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_trans_xbc_lo(const int i, const int j, const int k, const int n,
                                 const amrex::Array4<const amrex::Real> &s, 
                                       amrex::Real &lo,
                                       amrex::Real &hi,
                                       amrex::Real &uedge,
                                 const int bclo, 
                                 const int domlo[], const int domhi[], 
                                 bool ycouple, bool zcouple)
{
    using namespace amrex;

    // Low X
    if(ycouple && (j< domlo[1] || j > domhi[1])) return;
    if(zcouple && (k< domlo[2] || k > domhi[2])) return;

    if (i <= domlo[0])
    { 
        if(bclo==BCType::ext_dir){

#if 0
            if(n == XVEL){ //TODO not the right xvel? 
                if(uedge >= 0.e0){
                    lo = s(domlo[0]-1, j, k, n);
                    hi = s(domlo[0]-1, j, k, n);
                }
                else
                    lo = s(domlox, j, k, n);
            }
            else { 
#endif
                Real st = (uedge <= 1e-10) ? hi : s(domlo[0]-1,j,k,n);
                lo = st;
                hi = st;
//            }
        }
        else if(bclo == BCType::foextrap || bclo == BCType::hoextrap || bclo == BCType::reflect_even)
        {
            lo = hi;
        }
        else if(bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_trans_xbc_hi(const int i, const int j, const int k, const int n,
                                 const amrex::Array4<const amrex::Real> &s, 
                                       amrex::Real &lo,
                                       amrex::Real &hi,
                                       amrex::Real &uedge,
                                 const int bchi, 
                                 const int domlo[], const int domhi[], 
                                 bool ycouple, bool zcouple)
{
    using namespace amrex;

    if(ycouple && (j< domlo[1] || j > domhi[1])) return;
    if(zcouple && (k< domlo[2] || k > domhi[2])) return;

    // High X
    if (i >= domhi[0]) 
    {
        if(bchi==BCType::ext_dir){
#if 0
            if(n == XVEL){
                if(uedge <=0.){
                    lo = s(domhi[0],j,k,n);
                    hi = s(domhi[0],j,k,n);
                }  
                else  
                     hi = lo;
            }
            else { 
#endif
                amrex::Real st = (uedge >= -1e-10)? lo : s(domhi[0]+1,j,k,n);
                lo = st;
                hi = st;
//            }
        }
        else if (bchi == BCType::foextrap || bchi == BCType::hoextrap || bchi == BCType::reflect_even){
            hi = lo;
        }
        else if(bchi == BCType::reflect_odd){
            hi = 0.e0;
            lo = 0.e0;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_trans_ybc_lo(const int i, const int j, const int k, const int n,
                                 const amrex::Array4<const amrex::Real> &s, 
                                       amrex::Real &lo,
                                       amrex::Real &hi,
                                       amrex::Real vedge,
                                 const int bclo, 
                                 const int domlo[], const int domhi[], 
                                 bool xcouple, bool zcouple)
{
    using namespace amrex;

// Lower X 
    if(xcouple && (i < domlo[0] || i > domhi[0])) return;
    if(zcouple && (k < domlo[2] || k > domhi[2])) return;
    const int domloy = domlo[1];
    if(j <= domloy){ 
        if(bclo==BCType::ext_dir){
/*            if(n == YVEL){
                if(vedge >= 0.e0){
                    lo = s(i, domloy-1, k, n);
                    hi = s(i, domloy-1, k, n);
                } 
                else  
                    lo = s(i, domloy, k, n);
            }
            else{ */ 
                amrex::Real st = (vedge < 1e-10) ? hi : s(i,domloy-1,k,n);
                lo = st;
                hi = st;
//            }
        }
        else if(bclo == BCType::foextrap || bclo == BCType::hoextrap || bclo == BCType::reflect_even)
        {
            lo = hi;
        }
        else if(bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_trans_ybc_hi(const int i, const int j, const int k, const int n,
                                 const amrex::Array4<const amrex::Real> &s, 
                                       amrex::Real &lo,
                                       amrex::Real &hi,
                                       amrex::Real vedge,
                                 const int bchi, 
                                 const int domlo[], const int domhi[], 
                                 bool xcouple, bool zcouple)
{
    using namespace amrex;
    if(xcouple && (i < domlo[0] || i > domhi[0])) return;
    if(zcouple && (k < domlo[2] || k > domhi[2])) return;
    const int domhiy = domhi[1];
//Upper Y
    if(j > domhiy){
        if(bchi==BCType::ext_dir){
/*            if(n == YVEL){
                if(vedge <=0.){
                    lo = s(i,domhiy,k,n);
                    hi = s(i,domhiy,k,n);
                }
                else  
                hi = lo;
            } */ 
//            else{ 
                amrex::Real st = (vedge >= -1e-10)? lo : s(i,domhiy+1,k,n);
                lo = st;
                hi = st;
//            }
        }
        else if (bchi == BCType::foextrap || bchi == BCType::hoextrap || bchi == BCType::reflect_even){
            hi = lo;
        }
        else if(bchi == BCType::reflect_odd){
            hi = 0.e0;
            lo = 0.e0;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_trans_zbc_lo(const int i, const int j, const int k, const int n,
                                 const amrex::Array4<const amrex::Real> &s, 
                                       amrex::Real &lo,
                                       amrex::Real &hi,
                                       amrex::Real wedge,
                                 const int bclo, 
                                 const int domlo[], const int domhi[], 
                                 bool xcouple, bool ycouple)
{
    using namespace amrex;

    if(xcouple && (i < domlo[0] || i > domhi[0])) return;
    if(ycouple && (j < domlo[1] || j > domhi[1])) return;
    const int domloz = domlo[2];
// Lower Z
    if(k <= domloz){ 
        if(bclo==BCType::ext_dir){
/*            if(n == ZVEL){
                if(wedge >= 0.e0){
                    lo = s(i, j, domloz-1, n);
                    hi = s(i, j, domloz-1, n);
                }
                else  
                    lo = s(i, j, domloz, n);
            } 
            else{ */
                amrex::Real st = (wedge <= 1e-10) ? hi : s(i,j,domloz-1,n);
                lo = st;
                hi = st;
//            }
        }
        else if(bclo == BCType::foextrap || bclo == BCType::hoextrap || bclo == BCType::reflect_even)
        {
            lo = hi;
        }
        else if(bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_trans_zbc_hi(const int i, const int j, const int k, const int n,
                                 const amrex::Array4<const amrex::Real> &s, 
                                       amrex::Real &lo,
                                       amrex::Real &hi,
                                       amrex::Real wedge,
                                 const int bchi, 
                                 const int domlo[], const int domhi[], 
                                 bool xcouple, bool ycouple)
{
    using namespace amrex;
    if(xcouple && (i < domlo[0] || i > domhi[0])) return;
    if(ycouple && (j < domlo[1] || j > domhi[1])) return;
    const int domhiz = domhi[2];
//Upper Z 
    if(k > domhiz){
        if(bchi==BCType::ext_dir){
/*            if(n == ZVEL){
                if(wedge <=0.){
                    lo = s(i,j,domhiz+1,n);
                    hi = s(i,j,domhiz+1,n);
                }
                else  
                     hi = lo;
            }
            else{ */ 
                amrex::Real st = (wedge >= -1e-10)? lo : s(i,j,domhiz+1,n);
                lo = st;
                hi = st;
//            }
        }
        else if (bchi == BCType::foextrap || bchi == BCType::hoextrap || bchi == BCType::reflect_even){
            hi = lo;
        }
        else if(bchi == BCType::reflect_odd){
            hi = 0.e0;
            lo = 0.e0;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_cc_xbc(const int i, const int j, const int k, const int n,
                             const amrex::Array4<const amrex::Real> &s, 
                             amrex::Real &lo, amrex::Real &hi, 
                             const amrex::Array4< const amrex::Real> &uedge, 
                             const int bclo, const int bchi, 
                             const int domlo, const int domhi)
{
    using namespace amrex;
 
        if(i == domlo){
            if(bclo == BCType::ext_dir && uedge(domlo, j, k) >= 0.e0)
            {
                lo = s(domlo-1, j, k, n);
                hi = lo;
            }    
            else if(bclo == BCType::ext_dir && uedge(domlo, j, k) < 0)
                lo = hi;
            else if(bclo == BCType::foextrap || bclo == BCType::hoextrap)
                lo = hi;
            else if(bclo == BCType::reflect_even) lo = hi;
            else if(bclo == BCType::reflect_odd){
                hi = 0.e0;
                lo = 0.e0;
            }
        }
        else if(i==domhi+1){
            if(bchi == BCType::ext_dir && uedge(domhi+1, j, k) <= 0.e0)
            {
                lo = s(domhi+1, j, k, n);
                hi = lo;
            }    
            else if(bchi == BCType::ext_dir && uedge(domhi+1, j, k) > 0)
                hi = lo;
            else if(bchi == BCType::foextrap || bclo == BCType::hoextrap)
                hi = lo;
            else if(bchi == BCType::reflect_even) lo = hi;
            else if(bchi == BCType::reflect_odd){
                hi = 0.e0;
                lo = 0.e0;
            }
        }
        else return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_cc_ybc(const int i, const int j, const int k, const int n,
                             const amrex::Array4<const amrex::Real> &s, 
                             amrex::Real &lo, amrex::Real hi, 
                             const amrex::Array4< const amrex::Real> &vedge, 
                             const int bclo, const int bchi, 
                             const int domlo, const int domhi)
{
    using namespace amrex;

        if(j == domlo){
            if(bclo == BCType::ext_dir && vedge(i, domlo, k) >= 0.e0)
            {
                lo = s(i, domlo-1, k, n);
                hi = lo;
            }    
            else if(bclo == BCType::ext_dir && vedge(i, domlo, k) < 0)
                lo = hi;
            else if(bclo == BCType::foextrap || bclo == BCType::hoextrap)
                lo = hi;
            else if(bclo == BCType::reflect_even) lo = hi;
            else if(bclo == BCType::reflect_odd){
                hi = 0.e0;
                lo = 0.e0;
            }
        }
        else if(j==domhi+1){
            if(bchi == BCType::ext_dir && vedge(i,domhi+1, k) <= 0.e0)
            {
                lo = s(i, domhi+1, k, n);
                hi = lo;
            }    
            else if(bchi == BCType::ext_dir && vedge(i,domhi+1, k) > 0)
                hi = lo;
            else if(bchi == BCType::foextrap || bclo == BCType::hoextrap)
                hi = lo;
            else if(bchi == BCType::reflect_even) lo = hi;
            else if(bchi == BCType::reflect_odd){
                hi = 0.e0;
                lo = 0.e0;
            }
        }
        else return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_cc_zbc(const int i, const int j, const int k, const int n,
                             const amrex::Array4<const amrex::Real> &s, 
                             amrex::Real &lo, amrex::Real &hi, 
                             const amrex::Array4< const amrex::Real> &wedge, 
                             const int bclo, const int bchi, 
                             const int domlo, const int domhi)
{
    using namespace amrex;

        if(k == domlo){
            if(bclo == BCType::ext_dir && wedge(i,j, domlo) >= 0.e0)
            {
                lo = s(i, j, domlo-1, n);
                hi = lo;
            }    
            else if(bclo == BCType::ext_dir && wedge(i,j,domlo) < 0)
                lo = hi;
            else if(bclo == BCType::foextrap || bclo == BCType::hoextrap)
                lo = hi;
            else if(bclo == BCType::reflect_even) lo = hi;
            else if(bclo == BCType::reflect_odd){
                hi = 0.e0;
                lo = 0.e0;
            }
        }
        else if(k==domhi+1){
            if(bchi == BCType::ext_dir && wedge(i,j,domhi+1) <= 0.e0)
            {
                lo = s(i,j,domhi+1, n);
                hi = lo;
            }    
            else if(bchi == BCType::ext_dir && wedge(i,j, domhi+1) > 0)
                hi = lo;
            else if(bchi == BCType::foextrap || bclo == BCType::hoextrap)
                hi = lo;
            else if(bchi == BCType::reflect_even) lo = hi;
            else if(bchi == BCType::reflect_odd){
                hi = 0.e0;
                lo = 0.e0;
            }
        }
        else return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_corner_couple_yx (amrex::Real& lo1, amrex::Real& hi1,
                               int i, int j, int k, int n, amrex::Real dt, amrex::Real dx,
                               bool iconserv, amrex::Real lo, amrex::Real hi,
                               amrex::Array4<amrex::Real const> const& s,
                               amrex::Array4<amrex::Real const> const& divu_cc,
                               amrex::Array4<amrex::Real const> const& mac,
                               amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo - dt/(3.e0*dx)*(state(i+1,j-1,k,n)*mac(i+1,j-1,k)
                               - state(i  ,j-1,k,n)*mac(i  ,j-1,k))
            - dt/(3.e0)   *s(i,j-1,k,n)*divu_cc(i,j-1,k)
            + dt/(3.e0*dx)*s(i,j-1,k,n)*(mac(i+1,j-1,k)-mac(i,j-1,k));

        hi1  = hi - dt/(3.e0*dx)*(state(i+1,j,k,n)*mac(i+1,j,k)
                                - state(i  ,j,k,n)*mac(i  ,j,k))
            - dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
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
                               amrex::Array4<amrex::Real const> const& divu_cc,
                               amrex::Array4<amrex::Real const> const& mac,
                               amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo - dt/(3.e0*dx)*(state(i+1,j,k-1,n)*mac(i+1,j,k-1)
                               - state(i  ,j,k-1,n)*mac(i  ,j,k-1))
            - dt/(3.e0)   *s(i,j,k-1,n)*divu_cc(i,j,k-1)
            + dt/(3.e0*dx)*s(i,j,k-1,n)*(mac(i+1,j,k-1)-mac(i,j,k-1));

        hi1 = hi - dt/(3.e0*dx)*(state(i+1,j,k,n)*mac(i+1,j,k)
                               - state(i  ,j,k,n)*mac(i  ,j,k))
            - dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
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
                               amrex::Array4<amrex::Real const> const& divu_cc,
                               amrex::Array4<amrex::Real const> const& mac,
                               amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo - dt/(3.e0*dy)*(state(i-1,j+1,k,n)*mac(i-1,j+1,k)
                               - state(i-1,j  ,k,n)*mac(i-1,j  ,k))
            - dt/3.e0     *s(i-1,j,k,n)*divu_cc(i-1,j,k)
            + dt/(3.e0*dy)*s(i-1,j,k,n)*(mac(i-1,j+1,k)-mac(i-1,j,k));

        hi1 = hi - dt/(3.e0*dy)*(state(i,j+1,k,n)*mac(i,j+1,k)
                               - state(i,j  ,k,n)*mac(i,j  ,k))
            - dt/3.e0     *s(i,j,k,n)*divu_cc(i,j,k)
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
                               amrex::Array4<amrex::Real const> const& divu_cc,
                               amrex::Array4<amrex::Real const> const& mac,
                               amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo - dt/(3.e0*dy)*(state(i,j+1,k-1,n)*mac(i,j+1,k-1)
                               - state(i,j  ,k-1,n)*mac(i,j  ,k-1))
            - dt/(3.e0)   *s(i,j,k-1,n)*divu_cc(i,j,k-1)
            + dt/(3.e0*dy)*s(i,j,k-1,n)*(mac(i,j+1,k-1)-mac(i,j,k-1));

        hi1 = hi - dt/(3.e0*dy)*(state(i,j+1,k,n)*mac(i,j+1,k)
                               - state(i,j  ,k,n)*mac(i,j  ,k))
            - dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
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
                               amrex::Array4<amrex::Real const> const& divu_cc,
                               amrex::Array4<amrex::Real const> const& mac,
                               amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo - dt/(3.e0*dz)*(state(i-1,j,k+1,n)*mac(i-1,j,k+1)
                               - state(i-1,j,k  ,n)*mac(i-1,j,k))
            - dt/(3.e0)   *s(i-1,j,k,n)*divu_cc(i-1,j,k)
            + dt/(3.e0*dz)*s(i-1,j,k,n)*(mac(i-1,j,k+1) - mac(i-1,j,k));

        hi1 = hi - dt/(3.e0*dz)*(state(i,j,k+1,n)*mac(i,j,k+1)
                               - state(i,j,k  ,n)*mac(i,j,k))
            - dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
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
                               amrex::Array4<amrex::Real const> const& divu_cc,
                               amrex::Array4<amrex::Real const> const& mac,
                               amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo - dt/(3.e0*dz)*(state(i,j-1,k+1,n)*mac(i,j-1,k+1)
                               - state(i,j-1,k  ,n)*mac(i,j-1,k))
            - dt/(3.e0)   *s(i,j-1,k,n)*divu_cc(i,j-1,k)
            + dt/(3.e0*dz)*s(i,j-1,k,n)*(mac(i,j-1,k+1)-mac(i,j-1,k));

        hi1 = hi - dt/(3.e0*dz)*(state(i,j,k+1,n)*mac(i,j,k+1)
                               - state(i,j,k  ,n)*mac(i,j,k))
            - dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
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
