#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include "AMReX_RealVect.H"

#ifdef __cplusplus
extern "C"
{
#endif
    void set_delp_dir(int* delp_dir);

    void incflo_get_data(amrex::Real* gravity, amrex::Real* ro_0, amrex::Real* mu_0); 

    void incflo_set_cyclic(int* cyc_x, int* cyc_y, int* cyc_z);

    void set_bc_flow(amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength,
         amrex::Real* dx,      amrex::Real* dy,      amrex::Real* dz  );

    void set_bc_type(
        const int* bc_ilo, const int* bc_ihi,
	    const int* bc_jlo, const int* bc_jhi,
	    const int* bc_klo, const int* bc_khi,
        const int* domlo , const int* domhi,
        amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
        amrex::Real* xlen, amrex::Real* ylen, amrex::Real* zlen,
        const int* ng
    );

    void set_bc0(
        const int* slo, const int* shi,
        amrex::Real* ro,
        amrex::Real* mu,  amrex::Real* lambda,
        const int* bc_ilo, const int* bc_ihi,
	    const int* bc_jlo, const int* bc_jhi,
	    const int* bc_klo, const int* bc_khi,
        const int* domlo , const int* domhi,
        const int* ng    , const int* nodal_pressure
    );

    void init_fluid(const int* slo, const int* shi,
        const int*  lo, const int*  hi,
        const int* domlo , const int* domhi,
        amrex::Real* ro, amrex::Real* p,
        amrex::Real* vel,
        amrex::Real* mu,  amrex::Real* lambda,
        amrex::Real* dx,    amrex::Real* dy,       amrex::Real* dz,
        amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength
       );

    void init_fluid_restart( 
          const int* slo, const int* shi,
          const int*  lo, const int*  hi,
          amrex::Real* mu,    amrex::Real* lambda
       );

    void set_p0(const int*  lo   , const int*  hi,
	        const int* domlo , const int* domhi,
            amrex::Real* p0,  const int* rlo, const int* rhi,
            amrex::Real* gp0, const int* glo, const int* ghi,
	        amrex::Real* dx,       amrex::Real* dy,       amrex::Real* dz,
	        amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength,
	        int* delp_dir,
	        const int* bc_ilo, const int* bc_ihi,
        	const int* bc_jlo, const int* bc_jhi,
        	const int* bc_klo, const int* bc_khi,
            const int* ng, const int* nodal_pressure
       );

    void set_velocity_bcs (
        amrex::Real* vel, const int* ulo, const int* uhi,
	const int* bc_ilo, const int* bc_ihi,
	const int* bc_jlo, const int* bc_jhi,
	const int* bc_klo, const int* bc_khi,
	const int* domlo , const int* domhi,
        const int* ng, const int* extrap_dir_bcs
       );

    void set_vec_bcs (
        amrex::Real* vec, const int* ulo, const int* uhi,
	const int* bc_ilo, const int* bc_ihi,
	const int* bc_jlo, const int* bc_jhi,
	const int* bc_klo, const int* bc_khi,
	const int* domlo , const int* domhi,
        const int* ng
       );

    void zero_wall_norm_vel(
        const int* slo, const int* shi,
        amrex::Real* vel,
	const int* bc_ilo, const int* bc_ihi,
	const int* bc_jlo, const int* bc_jhi,
	const int* bc_klo, const int* bc_khi,
        const int* domlo , const int* domhi,
        const int* ng
        );

    void fill_bc0(amrex::Real* s, const int* slo, const int* shi,
	const int* bc_ilo, const int* bc_ihi,
	const int* bc_jlo, const int* bc_jhi,
	const int* bc_klo, const int* bc_khi,
        const int* domlo , const int* domhi,
        const int* ng
        );

    void init_periodic_vortices (
	const int* lo,  const int*  hi,
	amrex::Real* vel,  const int* ulo, const int* uhi,
	const amrex::Real* dx,
	const int*  domlo, const int* plane
       );

   
#ifdef __cplusplus
}
#endif
