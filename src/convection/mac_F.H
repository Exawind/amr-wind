#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void set_mac_velocity_bcs(amrex::Real* time, 
                              const int* slo, const int* shi,
                              amrex::Real* u,  const int* ulo, const int* uhi,
                              amrex::Real* v,  const int* vlo, const int* vhi,
                              amrex::Real* w,  const int* wlo, const int* whi,
                              const int* bc_ilo, const int* bc_ihi,
                              const int* bc_jlo, const int* bc_jhi,
                              const int* bc_klo, const int* bc_khi,
                              const int* domlo , const int* domhi,
                              const int* ng, const int* probtype);

#ifdef __cplusplus
			       }
#endif
