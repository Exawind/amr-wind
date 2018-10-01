#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void compute_viscosity(const int* lo, const int* hi,
                           amrex::Real* mu, const int* mulo, const int* muhi,
                           amrex::Real* strainrate, const int* slo, const int* shi,
                           const amrex::Real* dx);
#ifdef __cplusplus
}
#endif
