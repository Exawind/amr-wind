#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

   void add_forcing(
      const int* lo, const int* hi,
      amrex::Real* vel,  const int* ulo, const int* uhi,
      const amrex::Real* dt);

#ifdef __cplusplus
}
#endif
