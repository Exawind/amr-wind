#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void compute_mac_diveu (
	const int* lo, const int* hi,
	amrex::Real* diveu,  const int* slo, const int* shi,
	const amrex::Real* u, const int* ulo, const int* uhi,
	const amrex::Real* v, const int* vlo, const int* vhi,
	const amrex::Real* w, const int* wlo, const int* whi,
	const amrex::Real* dx
        );

   void compute_mac_diveu_eb (
      const int* lo, const int* hi,
      amrex::Real* diveu,  const int* slo, const int* shi,
      const amrex::Real* u, const int* ulo, const int* uhi,
      const amrex::Real* v, const int* vlo, const int* vhi,
      const amrex::Real* w, const int* wlo, const int* whi,
      const amrex::Real* afracx, const int* axlo, const int* axhi,
      const amrex::Real* afracy, const int* aylo, const int* ayhi,
      const amrex::Real* afracz, const int* azlo, const int* azhi,
      const amrex::Real* vfrac,  const int* vflo, const int* vfhi,
      const void*        flags,  const int*  flo, const int*  fhi, 
      const amrex::Real* dx
        );
   
    void project_mac_velocity (
	const int* lo, const int* hi,
	amrex::Real* u_i, const int* ulo, const int* uhi,
	const amrex::Real* dphidxi, const int* glo, const int* ghi,
	const amrex::Real* ro, const int* slo, const int* shi,
	const amrex::Real* c, const int* dir
        );

    void set_mac_velocity_bcs( 
        const int* slo, const int* shi,
        amrex::Real* u,  const int* ulo, const int* uhi,
        amrex::Real* v,  const int* vlo, const int* vhi,
        amrex::Real* w,  const int* wlo, const int* whi,
        const int* bc_ilo, const int* bc_ihi,
        const int* bc_jlo, const int* bc_jhi,
        const int* bc_klo, const int* bc_khi,
        const int* domlo , const int* domhi,
        const int* ng
        );

    void compute_bcoeff_mac (
	const int* lo, const int* hi,
	amrex::Real* bcoeff, const int* alo, const int* ahi,
	const amrex::Real* u,     const int* ulo, const int* uhi,
	const amrex::Real* ro, const int* slo, const int* shi,
	const int* dir
        );

#ifdef __cplusplus
			       }
#endif
