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
      const amrex::Real* etamax,
      const amrex::Real* gradp0max,
      const amrex::Real* dx,
      const amrex::Real* cfl,
      const         int* steady_state,
      const amrex::Real* time,
      const amrex::Real* stop_time,
      amrex::Real*  dt);
 
   void add_forcing (
      const int* lo, const int* hi,
      amrex::Real* vel,  const int* ulo, const int* uhi,
      const amrex::Real* dt);

#ifdef __cplusplus
}
#endif
