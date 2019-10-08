#ifndef __BC_MOD_F_H_
#define __BC_MOD_F_H_

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

int get_bc_defined(const int& pID);
void get_bc_type(const int& pID, char c_string[]);

amrex::Real get_bc_u(const int& pID);
amrex::Real get_bc_v(const int& pID);
amrex::Real get_bc_w(const int& pID);
amrex::Real get_bc_r(const int& pID);
amrex::Real get_bc_t(const int& pID);
amrex::Real get_bc_p(const int& pID);

amrex::Real get_delp_x();

amrex::Real get_delp_y();

amrex::Real get_delp_z();

void set_delp_x(const amrex::Real& delp_x_in);
void set_delp_y(const amrex::Real& delp_y_in);
void set_delp_z(const amrex::Real& delp_z_in);

int get_minf();
int get_pinf();
int get_pout();

void get_domain_bc (int* domain_bc_out);

#ifdef __cplusplus
}
#endif

#endif
