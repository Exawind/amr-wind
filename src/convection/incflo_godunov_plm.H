#ifndef GODUNOV_PLM_H
#define GODUNOV_PLM_H

#include <incflo_convection_K.H>
#include <AMReX_Gpu.H>
#include <AMReX_BCRec.H>

namespace {

// This version is called after the MAC projection
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Godunov_plm_fpu_x (const int i, const int j, const int k, const int n,
                        const amrex::Real dt, const amrex::Real dx,
                        amrex::Real& Im, amrex::Real& Ip,
                        const amrex::Array4<const amrex::Real> &S,
                        const amrex::Real& umac,
                        const amrex::BCRec bc,
                        const int domain_ilo, const int domain_ihi)
{
    using namespace amrex;
    {
        bool extdir_ilo = bc.lo(0) == BCType::ext_dir;
        bool extdir_ihi = bc.hi(0) == BCType::ext_dir;

        Real upls, umns;

        if (extdir_ilo && i == domain_ilo)
        {
            upls = S(i-1,j,k,n);
            umns = S(i-1,j,k,n);
        } else if (extdir_ihi && i == domain_ihi+1)
        {
            upls = S(i  ,j,k,n);
            umns = S(i  ,j,k,n);
        } else {
            upls = S(i  ,j,k,n) + 0.5 * (-1.0 - umac * dt/dx) * 
               incflo_ho_xslope_extdir(i  ,j,k,n,S, extdir_ilo, extdir_ihi, domain_ilo, domain_ihi);
            umns = S(i-1,j,k,n) + 0.5 * ( 1.0 - umac * dt/dx) * 
               incflo_ho_xslope_extdir(i-1,j,k,n,S, extdir_ilo, extdir_ihi, domain_ilo, domain_ihi);
        }

        Ip = umns;
        Im = upls;
    }
}

// This version is called after the MAC projection
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Godunov_plm_fpu_y (const int i, const int j, const int k, const int n,
                        const amrex::Real dt, const amrex::Real dy,
                        amrex::Real& Im, amrex::Real& Ip,
                        const amrex::Array4<const amrex::Real> &S,
                        const amrex::Real& vmac,
                        const amrex::BCRec bc,
                        const int domain_jlo, const int domain_jhi)
{
    using namespace amrex;
    {
        bool extdir_jlo = bc.lo(1) == BCType::ext_dir;
        bool extdir_jhi = bc.hi(1) == BCType::ext_dir;

        Real vpls, vmns;

        if (extdir_jlo && j == domain_jlo)
        {
            vpls = S(i,j-1,k,n);
            vmns = S(i,j-1,k,n);
        } else if (extdir_jhi && j == domain_jhi+1) {
            vpls = S(i,j  ,k,n);
            vmns = S(i,j  ,k,n);
        } else {
            vpls = S(i,j  ,k,n) + 0.5 * (-1.0 - vmac * dt/dy) * 
               incflo_ho_yslope_extdir(i,j  ,k,n,S, extdir_jlo, extdir_jhi, domain_jlo, domain_jhi);
            vmns = S(i,j-1,k,n) + 0.5 * ( 1.0 - vmac * dt/dy) * 
               incflo_ho_yslope_extdir(i,j-1,k,n,S, extdir_jlo, extdir_jhi, domain_jlo, domain_jhi);
        }

        Ip = vmns;
        Im = vpls;
    }
}

// This version is called after the MAC projection
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Godunov_plm_fpu_z (const int i, const int j, const int k, const int n,
                        const amrex::Real dt, const amrex::Real dz,
                        amrex::Real& Im, amrex::Real& Ip,
                        const amrex::Array4<const amrex::Real> &S,
                        const amrex::Real& wmac,
                        const amrex::BCRec bc,
                        const int domain_klo, const int domain_khi)
{
    using namespace amrex;
    {
        bool extdir_klo = bc.lo(2) == BCType::ext_dir;
        bool extdir_khi = bc.hi(2) == BCType::ext_dir;

        Real wpls, wmns;

        if (extdir_klo && k == domain_klo)
        {
            wpls = S(i,j,k-1,n);
            wmns = S(i,j,k-1,n);
        } else if (extdir_khi && k == domain_khi+1) {
            wpls = S(i,j,k  ,n);
            wmns = S(i,j,k  ,n);
        } else {
            wpls = S(i,j,k  ,n) + 0.5 * (-1.0 - wmac * dt/dz) * 
               incflo_ho_zslope_extdir(i,j,k  ,n,S, extdir_klo, extdir_khi, domain_klo, domain_khi);
            wmns = S(i,j,k-1,n) + 0.5 * ( 1.0 - wmac * dt/dz) * 
               incflo_ho_zslope_extdir(i,j,k-1,n,S, extdir_klo, extdir_khi, domain_klo, domain_khi);
        }

        Ip = wmns;
        Im = wpls;
    }
}


}
#endif 
