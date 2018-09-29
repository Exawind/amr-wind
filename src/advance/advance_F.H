#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

   void compute_new_dt ( 
      const amrex::Real* umax,
      const amrex::Real* vmax,
      const amrex::Real* wmax,
      const amrex::Real* romin,
      const amrex::Real* mumax,
      const amrex::Real* gradp0max,
      const amrex::Real* dx,
      const amrex::Real* cfl,
      const         int* steady_state,
      const amrex::Real* time,
      const amrex::Real* stop_time,
      amrex::Real*  dt);
 
   void compute_gradp0_max (
      const int* lo, const int* hi,
      amrex::Real* p0, const int* slo, const int* shi,
      const amrex::Real* gradp0max, const amrex::Real* dx,
      const int* nodal_pressure);

   void add_forcing (
      const int* lo, const int* hi,
      amrex::Real* vel,  const int* ulo, const int* uhi,
      const amrex::Real* ro, const int* slo, const int* shi,
      const int* domlo, const int* domhi,
      const amrex::Real* dx, const amrex::Real* dt);

   void compute_diveucc (
      const int* lo,  const int*  hi,
      amrex::Real* diveu,  const int* slo, const int* shi,
      const amrex::Real* vel,    const int* ulo, const int* uhi,
      const amrex::Real* dx);

   void compute_diveund (
      const int* lo,  const int*  hi,
      amrex::Real* diveu,  const int* slo, const int* shi,
      const amrex::Real* vec,    const int* ulo, const int* uhi,
      const amrex::Real* dx);
    
   void average_cc_to_fc ( 
      const int* lo, const int* hi,
      amrex::Real* fxc, const int* fxlo, const int* fxhi,
      amrex::Real* fyc, const int* fylo, const int* fyhi,
      amrex::Real* fzc, const int* fzlo, const int* fzhi,
      const amrex::Real* cc, const int* slo, const int* shi);

#ifdef __cplusplus
}
#endif
