#ifndef __PARAM_MOD_F_H_
#define __PARAM_MOD_F_H_

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

int get_dim_ic();

int get_dim_m();

int get_dim_bc();

amrex::Real get_my_huge();

int is_defined_db_cpp(const amrex::Real& x);

int is_undefined_db_cpp(const amrex::Real& x);

amrex::Real get_undefined();

#ifdef __cplusplus
}
#endif

AMREX_GPU_HOST_DEVICE
int is_equal(const amrex::Real& x, const amrex::Real& y);

#endif
