#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void compute_ugradu (
	const int* lo, const int* hi,
	amrex::Real* rhs, const int* rlo, const int* rhi, 
	const amrex::Real* vel   ,  const int* vello, const int* velhi,
	amrex::Real* u, const int* ulo, const int* uhi,
	amrex::Real* v, const int* vlo, const int* vhi,	
	amrex::Real* w, const int* wlo, const int* whi,
	const amrex::Real* xslopes, 
	const amrex::Real* yslopes, 
	const amrex::Real* zslopes, const int* slo, const int* shi, 
	const int* domlo, const int* domhi,
	const int* bc_ilo_type, const int* bc_ihi_type, 
	const int* bc_jlo_type, const int* bc_jhi_type,
	const int* bc_klo_type, const int* bc_khi_type,
	const amrex::Real* dx, const int* ng
        );


   void compute_ugradu_eb (
	const int* lo, const int* hi,
	amrex::Real* ugradu, const int* glo, const int* ghi, 
	const amrex::Real* vel   ,  const int* vello, const int* velhi,
	amrex::Real* u, const int* ulo, const int* uhi,
	amrex::Real* v, const int* vlo, const int* vhi,	
	amrex::Real* w, const int* wlo, const int* whi,
        const amrex::Real* afrac_x,const int*  axlo,const int*  axhi, 
        const amrex::Real* afrac_y,const int*  aylo,const int*  ayhi, 
        const amrex::Real* afrac_z,const int*  azlo,const int*  azhi, 
        const amrex::Real* cent_x,const int*   cxlo,const int*  cxhi, 
        const amrex::Real* cent_y,const int*   cylo,const int*  cyhi, 
        const amrex::Real* cent_z,const int*   czlo,const int*  czhi, 
        const void* flags, const int*  flo, const int*   fhi, 
        const amrex::Real* vfrac, const int* vflo, const int*  vfhi, 
        const amrex::Real* bcent, const int* blo, const int*  bhi, 
	const amrex::Real* xslopes, 
	const amrex::Real* yslopes, 
	const amrex::Real* zslopes, const int* slo, const int* shi, 
	const int* domlo, const int* domhi,
	const int* bc_ilo_type, const int* bc_ihi_type, 
	const int* bc_jlo_type, const int* bc_jhi_type,
	const int* bc_klo_type, const int* bc_khi_type,
	const amrex::Real* dx, const int* ng
      );

#ifdef __cplusplus
}
#endif
