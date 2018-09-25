#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include "AMReX_RealVect.H"

#ifdef BL_SINGLE_PRECISION_PARTICLES
typedef float amrex_particle_real;
#else
typedef double amrex_particle_real;
#endif

#ifdef __cplusplus
extern "C"
{
#endif

   void compute_new_dt ( 
      const amrex::Real* umax,
      const amrex::Real* vmax,
      const amrex::Real* wmax,
      const amrex::Real* romin,
      const amrex::Real* mumax,
      const amrex::Real* gradp0max,
      const amrex::Real* dx,
      const amrex::Real* cfl,
      const         int* steady_state,
      const amrex::Real* time,
      const amrex::Real* stop_time,
      amrex::Real*  dt);
 
   void compute_gradp0_max (
      const int* lo, const int* hi,
      amrex::Real* p0, const int* slo, const int* shi,
      const amrex::Real* gradp0max, const amrex::Real* dx,
      const int* nodal_pressure);

   void add_grad_pcc (
      const int* lo, const int* hi,
      amrex::Real* vel, const int* ulo, const int* uhi,
      const amrex::Real* ro, const int* slo, const int* shi,
      const amrex::Real* p, const amrex::Real* p0,
      const amrex::Real* dx, const amrex::Real* c);

   void add_grad_pnd (
      const int* lo, const int* hi,
      amrex::Real* vel, const int* ulo, const int* uhi,
      const amrex::Real* ro, const int* slo, const int* shi,
      const amrex::Real* p, 
      const amrex::Real* p0, const int* rlo, const int* rhi,
      const amrex::Real* dx, const amrex::Real* c);

   void add_grad_phicc (
      const int* lo, const int* hi,
      amrex::Real* vel, const int* ulo, const int* uhi,
      const amrex::Real* ro, const int* slo, const int* shi,
      const amrex::Real* phi, 
      const amrex::Real* dx, const amrex::Real* c);

   void add_grad_phind (
      const int* lo, const int* hi,
      amrex::Real* vel, const int* ulo, const int* uhi,
      const amrex::Real* ro, const int* slo, const int* shi,
      const amrex::Real* phi, const int* rlo, const int* rhi,
      const amrex::Real* dx, const amrex::Real* c);

   void compute_velocity_at_faces (
      const int* lo, const int* hi,
      const amrex::Real* u_mac, const int* ulo, const int* uhi, 
      const amrex::Real* v_mac, const int* vlo, const int* vhi, 
      const amrex::Real* w_mac, const int* wlo, const int* whi, 
      const amrex::Real* vel   ,  const int* vello, const int* velhi, 
      const amrex::Real* xslopes, const int* slo, const int* shi, 
      const amrex::Real* yslopes, 
      const amrex::Real* zslopes, 
      const int* bc_ilo_type, const int* bc_ihi_type, 
      const int* bc_jlo_type, const int* bc_jhi_type,
      const int* bc_klo_type, const int* bc_khi_type,
      const int* ng, 
      const int* domlo, const int* domhi );

   void compute_velocity_at_x_faces_eb (
      const int* lo, const int* hi,
            amrex::Real* u_mac,     const int* ulo,   const int* uhi, 
      const amrex::Real* vel   ,    const int* vello, const int* velhi, 
      const amrex::Real* slopes,    const int* slo,   const int* shi,
      const amrex::Real* areafrac,  const int* alo,   const int* ahi,
      const amrex::Real* cent,      const int* clo,   const int* chi,
      const void*        flags,     const int* flo,   const int* fhi,
      const int* bc_ilo_type,       const int* bc_ihi_type, 
      const int* ng, 
      const int* domlo, const int* domhi );

   void compute_velocity_at_y_faces_eb (
      const int* lo, const int* hi,
            amrex::Real* v_mac,     const int* vlo,   const int* vhi, 
      const amrex::Real* vel   ,    const int* vello, const int* velhi, 
      const amrex::Real* slopes,    const int* slo,   const int* shi,
      const amrex::Real* areafrac,  const int* alo,   const int* ahi,
      const amrex::Real* cent,      const int* clo,   const int* chi,
      const void*        flags,     const int* flo,   const int* fhi,
      const int* bc_jlo_type,       const int* bc_jhi_type, 
      const int* ng, 
      const int* domlo, const int* domhi );

   void compute_velocity_at_z_faces_eb (
      const int* lo, const int* hi,
            amrex::Real* w_mac,     const int* wlo,   const int* whi, 
      const amrex::Real* vel   ,    const int* vello, const int* velhi, 
      const amrex::Real* slopes,    const int* slo,   const int* shi,
      const amrex::Real* areafrac,  const int* alo,   const int* ahi,
      const amrex::Real* cent,      const int* clo,   const int* chi,
      const void*        flags,     const int* flo,   const int* fhi,
      const int* bc_klo_type,       const int* bc_khi_type, 
      const int* ng, 
      const int* domlo, const int* domhi );   

   void compute_slopes (
      const int* lo, const int* hi,
      const amrex::Real* u, const int* ulo, const int* uhi, 
      amrex::Real* xslopes, amrex::Real* yslopes,
      amrex::Real* zslopes, const int* slo, const int* shi, 
      const int* domlo, const int* domhi,
      const int* bc_ilo_type, const int* bc_ihi_type, 
      const int* bc_jlo_type, const int* bc_jhi_type,
      const int* bc_klo_type, const int* bc_khi_type,
      const int* nghost );

   void compute_slopes_eb (
      const int* lo, const int* hi,
      const amrex::Real* u, const int* ulo, const int* uhi, 
      amrex::Real* xslopes, amrex::Real* yslopes,
      amrex::Real* zslopes, const int* slo, const int* shi,
      const void* flags, const int* flo, const int* fhi,
      const int* domlo, const int* domhi,
      const int* bc_ilo_type, const int* bc_ihi_type, 
      const int* bc_jlo_type, const int* bc_jhi_type,
      const int* bc_klo_type, const int* bc_khi_type,
      const int* nghost );

   void set_ppe_bc (
      int* bc_lo, int* bc_hi,
      const int* domlo, const int* domhi,
      const int* ng, 
      const int* bct_ilo, const int* bct_ihi, 
      const int* bct_jlo, const int* bct_jhi,
      const int* bct_klo, const int* bct_khi);

   void compute_diveucc (
      const int* lo,  const int*  hi,
      amrex::Real* diveu,  const int* slo, const int* shi,
      const amrex::Real* vel,    const int* ulo, const int* uhi,
      const amrex::Real* dx);

   void compute_diveund (
      const int* lo,  const int*  hi,
      amrex::Real* diveu,  const int* slo, const int* shi,
      const amrex::Real* vec,    const int* ulo, const int* uhi,
      const amrex::Real* dx);
    
   void set_poisson_solver_bcs (
      const int* lo, const int* hi,
      amrex::Real* rhs, const int* slo, const int* shi,
      amrex::Real* phi, 
      const amrex::Real* oro_x, const int* ulo, const int* uhi,
      const amrex::Real* oro_y, const int* vlo, const int* vhi,
      const amrex::Real* oro_z, const int* wlo, const int* whi,
      const int* bc_lo, const int* bc_hi,
      const amrex::Real* dx, const int* domlo, const int* domhi );

   void add_forcing (
      const int* lo, const int* hi,
      amrex::Real* vel,  const int* ulo, const int* uhi,
      const amrex::Real* ro, const int* slo, const int* shi,
      const int* domlo, const int* domhi,
      const amrex::Real* dx, const amrex::Real* dt);
    
   void compute_bcoeff_cc (
      const int* lo,  const int*  hi,
      amrex::Real* bcoeff, const int* blo, const int* bhi,
      amrex::Real* ro, const int* slo, const int* shi,
      const int* dir );

   void compute_bcoeff_nd (
      const int* lo,  const int*  hi,
      amrex::Real* bcoeff, const int* blo, const int* bhi,
      const amrex::Real* ro, const int* slo, const int* shi,
      const int* dir );

   void set_phi (
      const int* lo, const int* hi,
      amrex::Real* phi, const int* slo, const int* shi,
      amrex::Real* pg,
      const amrex::Real* scale, const int* singular,
      const int*  domlo, const int *domhi );

   void extrap_pressure_to_ghost_cells (
      amrex::Real* p, const int* slo, const int* shi,
      const int* bc_ilo, const int* bc_ihi,
      const int* bc_jlo, const int* bc_jhi,
      const int* bc_klo, const int* bc_khi,
      const int* domlo , const int* domhi,
      const int* nghost);

   void set_scalar_bcs (
      amrex::Real* ro,
      amrex::Real* mu,    amrex::Real* lambda,
      const int* slo, const int* shi,
      const int* bc_ilo, const int* bc_ihi,
      const int* bc_jlo, const int* bc_jhi,
      const int* bc_klo, const int* bc_khi,
      const int* domlo , const int* domhi,
      const int* nghost);

   void set_pressure_bcs (
      amrex::Real* p, const int* slo, const int* shi,
      const int* bc_ilo, const int* bc_ihi,
      const int* bc_jlo, const int* bc_jhi,
      const int* bc_klo, const int* bc_khi,
      const int* domlo , const int* domhi);

   void average_cc_to_fc ( 
      const int* lo, const int* hi,
      amrex::Real* fxc, const int* fxlo, const int* fxhi,
      amrex::Real* fyc, const int* fylo, const int* fyhi,
      amrex::Real* fzc, const int* fzlo, const int* fzhi,
      const amrex::Real* cc, const int* slo, const int* shi);
#ifdef __cplusplus
}
#endif
