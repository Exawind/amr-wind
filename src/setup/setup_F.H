#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void set_delp_dir(int* delp_dir);
    
    void fortran_get_data
        (
            int* is_cyclic, 
            amrex::Real* delp,
            amrex::Real* gravity,
            amrex::Real* ro_0,
            amrex::Real* mu,
            amrex::Real* ic_u,
            amrex::Real* ic_v,
            amrex::Real* ic_w,
            amrex::Real* ic_p,
            amrex::Real* n,
            amrex::Real* tau_0,
            amrex::Real* papa_reg,
            amrex::Real* eta_0,
            const char* fluid_model_name, 
            int fluid_model_namelength
        ); 
    
    void set_bc_type
        (
            const int* bc_ilo, const int* bc_ihi,
            const int* bc_jlo, const int* bc_jhi,
            const int* bc_klo, const int* bc_khi,
            const int* domlo , const int* domhi,
            amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
            amrex::Real* xlen, amrex::Real* ylen, amrex::Real* zlen,
            const int* ng
        );
    
    void set_bc_mod
        (
            const int* pID, const int* pType,
            const amrex::Real* pLo, const amrex::Real* pHi,
            amrex::Real* pLoc,
            amrex::Real* pPg,
            amrex::Real* pVel
        );

    void init_fluid
        (
            const int* slo, const int* shi,
            const int*  lo, const int*  hi,
            const int* domlo, const int* domhi,
            amrex::Real* ro, 
            amrex::Real* p, 
            amrex::Real* vel,
            amrex::Real* eta, 
            amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
            amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength, 
            const int* probtype
        );
    
    void set_p0
        (
            const int* lo, const int* hi,
            const int* domlo, const int* domhi,
            amrex::Real* p0,  const int* rlo, const int* rhi,
            amrex::Real* gp0, 
            amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
            amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength,
            int* delp_dir,
            const int* bc_ilo, const int* bc_ihi,
        	const int* bc_jlo, const int* bc_jhi,
        	const int* bc_klo, const int* bc_khi,
            const int* ng
        );
   
#ifdef __cplusplus
}
#endif
