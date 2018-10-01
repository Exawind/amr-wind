#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

   void average_cc_to_fc ( 
      const int* lo, const int* hi,
      amrex::Real* fxc, const int* fxlo, const int* fxhi,
      amrex::Real* fyc, const int* fylo, const int* fyhi,
      amrex::Real* fzc, const int* fzlo, const int* fzhi,
      const amrex::Real* cc, const int* slo, const int* shi);

#ifdef __cplusplus
}
#endif
