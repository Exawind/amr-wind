#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void compute_vort (
	const int* lo,  const int*  hi,
	amrex::Real* vort , const int* slo, const int* shi,
	amrex::Real* vel, const int* ulo, const int* uhi,
	const amrex::Real* dx );
#ifdef __cplusplus
}
#endif
