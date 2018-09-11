#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include "AMReX_RealVect.H"

#ifdef BL_SINGLE_PRECISION_PARTICLES
typedef float amrex_particle_real;
#else
typedef double amrex_particle_real;
#endif

#ifdef __cplusplus
extern "C"
{
#endif

  void compute_divtau (
          const int* lo, const int* hi,
          amrex::Real* divtau, const int* dlo, const int* dhi,
          const amrex::Real* vel  , const int* vlo, const int* vhi,
          const amrex::Real* mu   ,  amrex::Real* lambda,
          const amrex::Real* ro   , const int* slo, const int* shi,
          const int* domlo, const int* domhi,
          const int* bc_ilo_type, const int* bc_ihi_type,
          const int* bc_jlo_type, const int* bc_jhi_type,
          const int* bc_klo_type, const int* bc_khi_type,
          const amrex::Real* dx, const int* ng,
          int* do_explicit_diffusion);

  void compute_divtau_eb (
          const int* lo, const int* hi,
          amrex::Real* divtau, const int* dlo, const int* dhi,
          const amrex::Real* vel  , const int* vlo, const int* vhi,
          const amrex::Real* mu   ,  amrex::Real* lambda,
          const amrex::Real* ro   , const int* slo, const int* shi,
          const void* flag, const int* fglo, const int* fghi,
          const int* domlo, const int* domhi,
          const int* bc_ilo_type, const int* bc_ihi_type,
          const int* bc_jlo_type, const int* bc_jhi_type,
          const int* bc_klo_type, const int* bc_khi_type,
          const amrex::Real* dx, const int* ng,
          int* do_explicit_diffusion);

  void set_diff_bc (
          int* bc_lo, int* bc_hi,
          const int* domlo, const int* domhi,
          const int* ng,
          const int* bct_ilo, const int* bct_ihi,
          const int* bct_jlo, const int* bct_jhi,
          const int* bct_klo, const int* bct_khi);

  void compute_bcoeff_diff (
          const int* lo,  const int*  hi,
          amrex::Real* bcoeff, const int* blo, const int* bhi,
          const amrex::Real* mu_g, const int* slo, const int* shi,
          const int* dir );

#ifdef __cplusplus
}
#endif
