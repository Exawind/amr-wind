#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void compute_divtau (
          const int* lo, const int* hi,
          amrex::Real* divtau, const int* dlo, const int* dhi,
          const amrex::Real* vel  , const int* vlo, const int* vhi,
          const amrex::Real* eta   ,
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
          const amrex::Real* eta   , 
          const amrex::Real* ro   , const int* slo, const int* shi,
          const void* flag, const int* fglo, const int* fghi,
          const amrex::Real* afrac_x,const int*  axlo,const int*  axhi, 
          const amrex::Real* afrac_y,const int*  aylo,const int*  ayhi, 
          const amrex::Real* afrac_z,const int*  azlo,const int*  azhi, 
          const amrex::Real* cent_x,const int*   cxlo,const int*  cxhi, 
          const amrex::Real* cent_y,const int*   cylo,const int*  cyhi, 
          const amrex::Real* cent_z,const int*   czlo,const int*  czhi,
          const amrex::Real* vfrac, const int* vflo, const int*  vfhi, 
          const amrex::Real* bcent, const int* blo, const int*  bhi, 
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

#ifdef __cplusplus
}
#endif
