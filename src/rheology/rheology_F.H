#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif
  AMREX_GPU_HOST_DEVICE
    amrex::Real viscosity(amrex::Real sr);
#ifdef __cplusplus
}
#endif
