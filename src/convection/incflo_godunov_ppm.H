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


AMREX_GPU_HOST_DEVICE 
inline amrex::Real vanLeer(const amrex::Real a, const amrex::Real b, const amrex::Real c)
{   
    amrex::Real dsc = 0.5e0*(b - c); 
    amrex::Real dsl = 2.e0*(a - c); 
    amrex::Real dsr = 2.e0*(b - a); 
    return (dsl*dsr > 0)? std::copysign(1, dsc)*amrex::min(std::abs(dsc), 
                                 amrex::min(std::abs(dsl), std::abs(dsr))) : 0.; 
}

AMREX_GPU_HOST_DEVICE
inline void Godunov_ppm_xbc(const int i, const int j, const int k, const int n, 
                           amrex::Real &sm, amrex::Real &sp, 
                           amrex::Real &sedge1, amrex::Real &sedge2, 
                           const amrex::Array4<const amrex::Real> &s, 
                           const int bclo, const int bchi, 
                           const int domlo, const int domhi)
{
        if(bclo == EXT_DIR || bclo == HOEXTRAP){
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

      if(bchi == EXT_DIR || bchi == HOEXTRAP){
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

AMREX_GPU_HOST_DEVICE
inline void Godunov_ppm_ybc(const int i, const int j, const int k, const int n, 
                           amrex::Real &sm, amrex::Real &sp, 
                           amrex::Real &sedge1, amrex::Real &sedge2, 
                           const amrex::Array4<const amrex::Real> &s, 
                           const int bclo, const int bchi, 
                           const int domlo, const int domhi)
{
        if(bclo == EXT_DIR || bclo == HOEXTRAP){
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

      if(bchi == EXT_DIR || bchi == HOEXTRAP){
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

AMREX_GPU_HOST_DEVICE
inline void Godunov_ppm_zbc(const int i, const int j, const int k, const int n, 
                           amrex::Real &sm, amrex::Real &sp, 
                           amrex::Real &sedge1, amrex::Real &sedge2, 
                           const amrex::Array4<const amrex::Real> &s, 
                           const int bclo, const int bchi, 
                           const int domlo, const int domhi)
{
        if(bclo == EXT_DIR || bclo == HOEXTRAP){
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

      if(bchi == EXT_DIR || bchi == HOEXTRAP){
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
AMREX_GPU_HOST_DEVICE
inline void Godunov_ppm(const int i, const int j, const int k, const int n, 
                        const amrex::Real dt, const amrex::Real dx, 
                        const amrex::Array4<const amrex::Real> &S, 
                        const amrex::Array4<const amrex::Real> &vel_edge, 
                        const amrex::Array4<amrex::Real> &Im, 
                        const amrex::Array4<amrex::Real> &Ip,
                        const amrex::BCRec bc, const int domlo, const int domhi,                  
                        const int dir)
{
    amrex::Real sedge1, sedge2, sp, sm; 
    amrex::Real s1 = S(i,j,k,n);     
    
    if (dir==0)
    {
        amrex::Real s2 = S(i-1,j,k,n); 
        amrex::Real d1 = vanLeer(S(i,j,k,n), S(i+1,j,k,n), S(i-1,j,k,n)); 
        amrex::Real d2 = vanLeer(S(i-1,j,k,n), S(i,j,k,n), S(i-2,j,k,n)); 
        sedge1 = 0.5e0*(S(i,j,k,n) + S(i-1,j,k,n)) - 1.e0/6.e0*(d1 - d2); 
        sedge1 = amrex::max(sedge1, amrex::min(S(i,j,k,n), S(i-1,j,k,n))); 
        sedge1 = amrex::min(sedge1, amrex::max(S(i,j,k,n), S(i-1,j,k,n))); 

        amrex::Real s3 = s2;
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
            Ip(i,j,k,n) = sp - (0.5*sigmap)*((sp - sm) 
                             - (1.e0 -2.e0/3.e0*sigmap)*s6);
        }
        else
            Ip(i,j,k,n) = S(i,j,k,n); 

        if(vel_edge(i,j,k) < -1.e-06)
            Im(i,j,k,n) = sm + (0.5*sigmam)*((sp-sm)
                             + (1.e0 - 2.e0/3.e0*sigmam)*s6); 
        else
            Im(i,j,k,n) = S(i,j,k,n); 
    }
    else if (dir==1)
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
            Ip(i,j,k,n) = sp - (0.5*sigmap)*((sp - sm) 
                             - (1.e0 -2.e0/3.e0*sigmap)*s6); 
        else
            Ip(i,j,k,n) = S(i,j,k,n); 

        if(vel_edge(i,j,k) < -1.e-06)
            Im(i,j,k,n) = sm + (0.5*sigmam)*((sp-sm) 
                             + (1.e0 - 2.e0/3.e0*sigmam)*s6); 
        else
            Im(i,j,k,n) = S(i,j,k,n); 

    }    
    else
    {
        amrex::Real s2 = S(i,j,k-1,n); 

        amrex::Real d1 = vanLeer(S(i,j,k,n), S(i,j,k+1,n), S(i,j,k-1,n)); 
        amrex::Real d2 = vanLeer(S(i,j,k-1,n), S(i,j,k,n), S(i,j,k-2,n)); 

        sedge1 = 0.5e0*(S(i,j,k,n) + S(i,j,k-1,n)) - 1.e0/6.e0*(d1 - d2); 
        sedge1 = amrex::max(sedge1, amrex::min(S(i,j,k,n), S(i,j,k-1,n))); 
        sedge1 = amrex::min(sedge1, amrex::max(S(i,j,k,n), S(i,j,k-1,n))); 

        amrex::Real s3 = s2;
        s2 = s1; 
        s1 = S(i,j,k+1,n); 

        amrex::Real d3 = d2;         
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
            Ip(i,j,k,n) = sp - (0.5*sigmap)*((sp - sm)
                             - (1.e0 -2.e0/3.e0*sigmap)*s6); 
        else
            Ip(i,j,k,n) = s2;


        if(vel_edge(i,j,k) < -1.e-06)
            Im(i,j,k,n) = sm + (0.5*sigmam)*((sp-sm) 
                             + (1.e0 - 2.e0/3.e0*sigmam)*s6); 
        else
            Im(i,j,k,n) = s2;
    }
}


AMREX_GPU_HOST_DEVICE
inline void Godunov_trans_xbc(const int i, const int j, const int k, const int n, 
                             const amrex::Array4<const amrex::Real> &s, 
                             amrex::Real &lo, amrex::Real &hi, 
                             const amrex::Array4< const amrex::Real> &uedge, 
                             const int bclo, const int bchi, 
                             const int domlo[], const int domhi[], 
                             bool ycouple, bool zcouple)
{
// Lower X
    if(ycouple && (j< domlo[1] || j > domhi[1])) return; 
    if(zcouple && (k< domlo[2] || k > domhi[2])) return;
    const int domlox = domlo[0], domhix = domhi[0]; 
    if(i <= domlox){ 
        if(bclo==EXT_DIR){
/*            if(n == XVEL){ //TODO not the right xvel? 
                if(uedge(domlox,j,k) >= 0.e0){
                    lo = s(domlox-1, j, k, n); 
                    hi = s(domlox-1, j, k, n); 
                }
                else
                    lo = s(domlox, j, k, n); 
            }
            else{ */
                amrex::Real st = (uedge(domlox,j,k) <= 1e-10) ? hi : s(domlox-1,j,k,n); 
                lo = st; 
                hi = st; 
//            }
        }
        else if(bclo == FOEXTRAP || bclo == HOEXTRAP || bclo == REFLECT_EVEN)
        {
            lo = hi; 
        }
        else if(bclo == REFLECT_ODD)
        {
            hi = 0.; 
            lo = 0.; 
        }
    }
//Upper X 
    else if(i >= domhix){
        if(bchi==EXT_DIR){
/*            if(n == XVEL){
                if(uedge(domhix,j,k) <=0.){
                    lo = s(domhix,j,k,n); 
                    hi = s(domhix,j,k,n);
                }  
                else  
                     hi = lo; 
            }
            else{ */ 
                amrex::Real st = (uedge(domhix+1,j,k) >= -1e-10)? lo : s(domhix+1,j,k,n); 
                lo = st; 
                hi = st; 
//            }
        }
        else if (bchi == FOEXTRAP || bchi == HOEXTRAP || bchi == REFLECT_EVEN){
            hi = lo; 
        }
        else if(bchi == REFLECT_ODD){
            hi = 0.e0; 
            lo = 0.e0; 
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline void Godunov_trans_ybc(const int i, const int j, const int k, const int n, 
                             const amrex::Array4<const amrex::Real> &s, 
                             amrex::Real &lo, amrex::Real hi, 
                             const amrex::Array4< const amrex::Real> &vedge, 
                             const int bclo, const int bchi, 
                             const int domlo[], const int domhi[], 
                             bool xcouple, bool zcouple)
{
// Lower X 
    if(xcouple && (i < domlo[0] || i > domhi[0])) return; 
    if(zcouple && (k < domlo[2] || k > domhi[2])) return; 
    const int domloy = domlo[1], domhiy = domhi[1]; 
    if(j <= domloy){ 
        if(bclo==EXT_DIR){
/*            if(n == YVEL){
                if(vedge(i,domloy,k) >= 0.e0){
                    lo = s(i, domloy-1, k, n); 
                    hi = s(i, domloy-1, k, n); 
                } 
                else  
                    lo = s(i, domloy, k, n); 
            }
            else{ */ 
                amrex::Real st = (vedge(i,domloy,k) <= 1e-10) ? hi : s(i,domloy-1,k,n); 
                lo = st; 
                hi = st; 
//            }
        }
        else if(bclo == FOEXTRAP || bclo == HOEXTRAP || bclo == REFLECT_EVEN)
        {
            lo = hi; 
        }
        else if(bclo == REFLECT_ODD)
        {
            hi = 0.; 
            lo = 0.; 
        }
    }
//Upper Y
    else if(j > domhiy){
        if(bchi==EXT_DIR){
/*            if(n == YVEL){
                if(vedge(i,domhiy,k) <=0.){
                    lo = s(i,domhiy,k,n); 
                    hi = s(i,domhiy,k,n);
                }
                else  
                hi = lo; 
            } */ 
//            else{ 
                amrex::Real st = (vedge(i,domhiy+1,k) >= -1e-10)? lo : s(i,domhiy+1,k,n); 
                lo = st; 
                hi = st; 
//            }
        }
        else if (bchi == FOEXTRAP || bchi == HOEXTRAP || bchi == REFLECT_EVEN){
            hi = lo; 
        }
        else if(bchi == REFLECT_ODD){
            hi = 0.e0; 
            lo = 0.e0; 
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline void Godunov_trans_zbc(const int i, const int j, const int k, const int n, 
                             const amrex::Array4<const amrex::Real> &s, 
                             amrex::Real &lo, amrex::Real &hi, 
                             const amrex::Array4< const amrex::Real> &wedge, 
                             const int bclo, const int bchi, 
                             const int domlo[], const int domhi[],
                             bool xcouple, bool ycouple)
{
    if(xcouple && (i < domlo[0] || i > domhi[0])) return; 
    if(ycouple && (j < domlo[1] || j > domhi[1])) return; 
    const int domloz = domlo[2], domhiz = domhi[2]; 
// Lower Z
    if(k <= domloz){ 
        if(bclo==EXT_DIR){
/*            if(n == ZVEL){
                if(wedge(i,j,domloz) >= 0.e0){
                    lo = s(i, j, domloz-1, n); 
                    hi = s(i, j, domloz-1, n); 
                }
                else  
                    lo = s(i, j, domloz, n); 
            } 
            else{ */
                amrex::Real st = (wedge(i,j,domloz) <= 1e-10) ? hi : s(i,j,domloz-1,n); 
                lo = st; 
                hi = st;
//            }
        }
        else if(bclo == FOEXTRAP || bclo == HOEXTRAP || bclo == REFLECT_EVEN)
        {
            lo = hi; 
        }
        else if(bclo == REFLECT_ODD)
        {
            hi = 0.; 
            lo = 0.; 
        }
    }
//Upper Z 
    else if(k > domhiz){
        if(bchi==EXT_DIR){
/*            if(n == ZVEL){
                if(wedge(i,j,domhiz) <=0.){
                    lo = s(i,j,domhiz+1,n); 
                    hi = s(i,j,domhiz+1,n);
                }
                else  
                     hi = lo; 
            }
            else{ */ 
                amrex::Real st = (wedge(i,j,domhiz+1) >= -1e-10)? lo : s(i,j,k,domhiz+1); 
                lo = st; 
                hi = st; 
//            }
        }
        else if (bchi == FOEXTRAP || bchi == HOEXTRAP || bchi == REFLECT_EVEN){
            hi = lo; 
        }
        else if(bchi == REFLECT_ODD){
            hi = 0.e0; 
            lo = 0.e0; 
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline void Godunov_cc_xbc(const int i, const int j, const int k, const int n, 
                             const amrex::Array4<const amrex::Real> &s, 
                             amrex::Real &lo, amrex::Real &hi, 
                             const amrex::Array4< const amrex::Real> &uedge, 
                             const int bclo, const int bchi, 
                             const int domlo, const int domhi)
{ 
        if(i == domlo){
            if(bclo == EXT_DIR && uedge(domlo, j, k) >= 0.e0)
            {
                lo = s(domlo-1, j, k, n); 
                hi = lo; 
            }    
            else if(bclo == EXT_DIR && uedge(domlo, j, k) < 0)
                lo = hi; 
            else if(bclo == FOEXTRAP || bclo == HOEXTRAP)
                lo = hi; 
            else if(bclo == REFLECT_EVEN) lo = hi; 
            else if(bclo == REFLECT_ODD){
                hi = 0.e0; 
                lo = 0.e0; 
            }
        }
        else if(i==domhi+1){
            if(bchi == EXT_DIR && uedge(domhi+1, j, k) <= 0.e0)
            {
                lo = s(domhi+1, j, k, n); 
                hi = lo; 
            }    
            else if(bchi == EXT_DIR && uedge(domhi+1, j, k) > 0)
                hi = lo; 
            else if(bchi == FOEXTRAP || bclo == HOEXTRAP)
                hi = lo; 
            else if(bchi == REFLECT_EVEN) lo = hi; 
            else if(bchi == REFLECT_ODD){
                hi = 0.e0; 
                lo = 0.e0; 
            }
        }
        else return; 
}

AMREX_GPU_HOST_DEVICE
inline void Godunov_cc_ybc(const int i, const int j, const int k, const int n, 
                             const amrex::Array4<const amrex::Real> &s, 
                             amrex::Real &lo, amrex::Real hi, 
                             const amrex::Array4< const amrex::Real> &vedge, 
                             const int bclo, const int bchi, 
                             const int domlo, const int domhi)
{
        if(j == domlo){
            if(bclo == EXT_DIR && vedge(i, domlo, k) >= 0.e0)
            {
                lo = s(i, domlo-1, k, n); 
                hi = lo; 
            }    
            else if(bclo == EXT_DIR && vedge(i, domlo, k) < 0)
                lo = hi; 
            else if(bclo == FOEXTRAP || bclo == HOEXTRAP)
                lo = hi; 
            else if(bclo == REFLECT_EVEN) lo = hi; 
            else if(bclo == REFLECT_ODD){
                hi = 0.e0; 
                lo = 0.e0; 
            }
        }
        else if(j==domhi+1){
            if(bchi == EXT_DIR && vedge(i,domhi+1, k) <= 0.e0)
            {
                lo = s(i, domhi+1, k, n); 
                hi = lo; 
            }    
            else if(bchi == EXT_DIR && vedge(i,domhi+1, k) > 0)
                hi = lo; 
            else if(bchi == FOEXTRAP || bclo == HOEXTRAP)
                hi = lo; 
            else if(bchi == REFLECT_EVEN) lo = hi; 
            else if(bchi == REFLECT_ODD){
                hi = 0.e0; 
                lo = 0.e0; 
            }
        }
        else return; 
}

AMREX_GPU_HOST_DEVICE
inline void Godunov_cc_zbc(const int i, const int j, const int k, const int n, 
                             const amrex::Array4<const amrex::Real> &s, 
                             amrex::Real &lo, amrex::Real &hi, 
                             const amrex::Array4< const amrex::Real> &wedge, 
                             const int bclo, const int bchi, 
                             const int domlo, const int domhi)
{
        if(k == domlo){
            if(bclo == EXT_DIR && wedge(i,j, domlo) >= 0.e0)
            {
                lo = s(i, j, domlo-1, n); 
                hi = lo; 
            }    
            else if(bclo == EXT_DIR && wedge(i,j,domlo) < 0)
                lo = hi; 
            else if(bclo == FOEXTRAP || bclo == HOEXTRAP)
                lo = hi; 
            else if(bclo == REFLECT_EVEN) lo = hi; 
            else if(bclo == REFLECT_ODD){
                hi = 0.e0; 
                lo = 0.e0; 
            }
        }
        else if(k==domhi+1){
            if(bchi == EXT_DIR && wedge(i,j,domhi+1) <= 0.e0)
            {
                lo = s(i,j,domhi+1, n); 
                hi = lo; 
            }    
            else if(bchi == EXT_DIR && wedge(i,j, domhi+1) > 0)
                hi = lo; 
            else if(bchi == FOEXTRAP || bclo == HOEXTRAP)
                hi = lo; 
            else if(bchi == REFLECT_EVEN) lo = hi; 
            else if(bchi == REFLECT_ODD){
                hi = 0.e0; 
                lo = 0.e0; 
            }
        }
        else return; }




AMREX_GPU_HOST_DEVICE
inline void Godunov_corner_couple(const int i, const int j, const int k, const int n, 
                                 const amrex::Real dt, const amrex::Real dx, 
                                       int iconserv[], 
                                 const amrex::Array4<const amrex::Real> &lo, 
                                 const amrex::Array4<const amrex::Real> &hi,  
                                 const amrex::Array4<const amrex::Real> &s, 
                                 const amrex::Array4<const amrex::Real> &divu, 
                                 const amrex::Array4<const amrex::Real> &mac,         
                                 const amrex::Array4<const amrex::Real> &state, 
                                 const amrex::Array4<amrex::Real> &lo1, 
                                 const amrex::Array4<amrex::Real> &hi1,
                                 const int dir1, const int dir2)
{
//=========================== X ====================================================
    if(dir1 ==0){
        if(iconserv[n]==1){
        //xylo xyhi 
            if(dir2 == 0){
                lo1(i,j,k,n)  = lo(i, j, k, n)
                              - dt/(3.e0*dx)*(state(i-1,j+1,k,n)*mac(i-1,j+1,k)
                              - state(i-1,j,k,n)*mac(i-1,j,k))
                              - dt/3.e0*s(i-1,j,k,n)*divu(i-1,j,k)
                              + dt/(3.e0*dx)*s(i-1,j,k,n)*(mac(i-1,j+1,k)-mac(i-1,j,k));

                hi1(i,j,k,n)  = hi(i, j, k, n)
                              - dt/(3.e0*dx)*(state(i,j+1,k,n)*mac(i,j+1,k)
                              - state(i,j,k,n)*mac(i,j,k))
                              - dt/3.e0*s(i,j,k,n)*divu(i,j,k)
                              + dt/(3.e0*dx)*s(i,j,k,n)*(mac(i,j+1,k)-mac(i,j,k));
            }
        //xzlo xzhi 
            if(dir2 == 1){
                lo1(i,j,k,n)  = lo(i,j,k,n)
                              - dt/(3.e0*dx)*(state(i-1,j,k+1,n)*mac(i-1,j,k+1)
                              - state(i-1,j,k,n)*mac(i-1,j,k))
                              - dt/(3.e0)*s(i-1,j,k,n)*divu(i-1,j,k)
                              + dt/(3.e0*dx)*s(i-1,j,k,n)*(mac(i-1,j,k+1) - mac(i-1,j,k));

                hi1(i,j,k,n)  = hi(i,j,k,n)
                              - dt/(3.e0*dx)*(state(i,j,k+1,n)*mac(i,j,k+1)
                              - state(i,j,k,n)*mac(i,j,k))
                              - dt/(3.e0)*s(i,j,k,n)*divu(i,j,k)
                              + dt/(3.e0*dx)*s(i,j,k,n)*(mac(i,j,k+1) - mac(i,j,k));
            }
            }
        else{
            if(dir2 == 0){
               lo1(i,j,k,n)  = lo(i,j,k,n)
                             - dt/(6.e0*dx)*(mac(i-1,j+1,k) + mac(i-1,j,k))
                             * (state(i-1,j+1,k,n) - state(i-1,j,k,n));
               hi1(i,j,k,n)  = hi(i,j,k,n)
                             - dt/(6.e0*dx)*(mac(i,j+1,k) + mac(i,j,k))
                             * (state(i, j+1,k, n) - state(i, j, k, n));
            }
            if(dir2 == 1){
               lo1(i,j,k,n)  = lo(i,j,k,n)
                             - dt/(6.e0*dx)*(mac(i-1,j,k+1) + mac(i-1,j,k))
                             * (state(i-1,j,k+1,n) - state(i-1,j,k,n)); 
                
               hi1(i,j,k,n)  = hi(i,j,k,n) 
                             - dt/(6.e0*dx)*(mac(i,j,k+1) + mac(i,j,k))
                             * (state(i,j,k+1,n) - state(i,j,k,n)); 
            }
        }
    }
//=========================== Y ========================================================
    if(dir1==1){
        if(iconserv[n]==1){
        //yxlo yxhi 
             if(dir2 == 0){
   
                lo1(i,j,k,n)  = lo(i,j,k,n) 
                              - dt/(3.e0*dx)*(state(i+1,j-1,k,n)*mac(i+1,j-1,k)
                              - state(i,j-1,k,n)*mac(i,j-1,k))
                              - dt/(3.e0)*s(i,j-1,k,n)*divu(i,j-1,k)
                              + dt/(3.e0*dx)*s(i,j-1,k,n)*(mac(i+1,j-1,k)-mac(i,j-1,k));

                hi1(i,j,k,n)  = hi(i,j,k,n) 
                              - dt/(3.e0*dx)*(state(i+1,j,k,n)*mac(i+1,j,k)
                              - state(i,j,k,n)*mac(i,j,k))
                              - dt/(3.e0)*s(i,j,k,n)*divu(i,j,k) 
                              + dt/(3.e0*dx)*s(i,j,k,n)*(mac(i+1,j,k) - mac(i,j,k));
            }
        //yzlo yzhi         
             if(dir2 == 1){
               lo1(i,j,k,n)  = lo(i,j,k,n) 
                              - dt/(3.e0*dx)*(state(i,j-1,k+1,n)*mac(i,j-1,k+1)
                              - state(i,j-1,k,n)*mac(i,j-1,k))
                              - dt/(3.e0)*s(i,j-1,k,n)*divu(i,j-1,k)
                              + dt/(3.e0*dx)*s(i,j-1,k,n)*(mac(i,j-1,k+1)-mac(i,j-1,k));

                hi1(i,j,k,n)  = hi(i,j,k,n) 
                              - dt/(3.e0*dx)*(state(i,j,k+1,n)*mac(i,j,k+1)
                              - state(i,j,k,n)*mac(i,j,k))
                              - dt/(3.e0)*s(i,j,k,n)*divu(i,j,k) 
                              + dt/(3.e0*dx)*s(i,j,k,n)*(mac(i,j,k+1) - mac(i,j,k));
            }
        }
        else{
             if(dir2 == 0){
              lo1(i,j,k,n)  = lo(i,j,k,n)
                             - dt/(6.e0*dx)*(mac(i+1,j-1,k) + mac(i,j-1,k))
                             * (state(i+1,j-1,k,n) - state(i,j-1,k,n));
               hi1(i,j,k,n)  = hi(i,j,k,n)
                             - dt/(6.e0*dx)*(mac(i+1,j,k) + mac(i,j,k))
                             * (state(i+1, j,k, n) - state(i, j, k, n));
            }
             if(dir2 == 1){
              lo1(i,j,k,n)  = lo(i,j,k,n)
                             - dt/(6.e0*dx)*(mac(i,j-1,k+1) + mac(i,j-1,k))
                             * (state(i,j-1,k+1,n) - state(i,j-1,k,n)); 
                
               hi1(i,j,k,n)  = hi(i,j,k,n) 
                             - dt/(6.e0*dx)*(mac(i,j,k+1) + mac(i,j,k))
                             * (state(i,j,k+1,n) - state(i,j,k,n)); 
            }
        }
    }
//=========================== Z ======================================================    
    if(dir1==2){
        if(iconserv[n]==1){
        //zxlo zxhi 
             if(dir2 == 0){
               lo1(i,j,k,n)  = lo(i,j,k,n) 
                              - dt/(3.e0*dx)*(state(i+1,j,k-1,n)*mac(i+1,j,k-1)
                              - state(i,j,k-1,n)*mac(i,j,k-1))
                              - dt/(3.e0)*s(i,j,k-1,n)*divu(i,j,k-1)
                              + dt/(3.e0*dx)*s(i,j,k-1,n)*(mac(i+1,j,k-1)-mac(i,j,k-1));

                hi1(i,j,k,n)  = hi(i,j,k,n) 
                              - dt/(3.e0*dx)*(state(i+1,j,k,n)*mac(i+1,j,k)
                              - state(i,j,k,n)*mac(i,j,k))
                              - dt/(3.e0)*s(i,j,k,n)*divu(i,j,k) 
                              + dt/(3.e0*dx)*s(i,j,k,n)*(mac(i+1,j,k) - mac(i,j,k));
            }
        //zylo zyhi         
             if(dir2 == 1){
               lo1(i,j,k,n)  = lo(i,j,k,n) 
                              - dt/(3.e0*dx)*(state(i,j+1,k-1,n)*mac(i,j+1,k-1)
                              - state(i,j,k-1,n)*mac(i,j,k-1))
                              - dt/(3.e0)*s(i,j,k-1,n)*divu(i,j,k-1)
                              + dt/(3.e0*dx)*s(i,j,k-1,n)*(mac(i,j+1,k-1)-mac(i,j,k-1));

                hi1(i,j,k,n)  = hi(i,j,k,n) 
                              - dt/(3.e0*dx)*(state(i,j+1,k,n)*mac(i,j+1,k)
                              - state(i,j,k,n)*mac(i,j,k))
                              - dt/(3.e0)*s(i,j,k,n)*divu(i,j,k) 
                              + dt/(3.e0*dx)*s(i,j,k,n)*(mac(i,j+1,k) - mac(i,j,k));
            }
        }
        else{
             if(dir2 == 0){
              lo1(i,j,k,n)  = lo(i,j,k,n)
                             - dt/(6.e0*dx)*(mac(i+1,j,k-1) + mac(i,j,k-1))
                             * (state(i+1,j,k-1,n) - state(i,j,k-1,n));
               hi1(i,j,k,n)  = hi(i,j,k,n)
                             - dt/(6.e0*dx)*(mac(i+1,j,k) + mac(i,j,k))
                             * (state(i+1, j,k, n) - state(i, j, k, n));
             }
             if(dir2 == 1){
              lo1(i,j,k,n)  = lo(i,j,k,n)
                             - dt/(6.e0*dx)*(mac(i,j+1,k-1) + mac(i,j,k-1))
                             * (state(i,j+1,k-1,n) - state(i,j,k-1,n)); 
                
               hi1(i,j,k,n)  = hi(i,j,k,n) 
                             - dt/(6.e0*dx)*(mac(i,j+1,k) + mac(i,j,k))
                             * (state(i,j+1,k,n) - state(i,j,k,n)); 
            }
        }
    }

}
#endif 
