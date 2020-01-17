#include "incflo_godunov_ppm.H" 
#include "incflo.H"

using namespace amrex;

#if 0
void
incflo::incflo_godunov_fluxes_on_box (const int lev, Box& bx,
                                      const Array4<Real> &a_fx,
                                      const Array4<Real> &a_fy, 
                                      const Array4<Real> &a_fz, 
                                      const Array4<Real> &s, 
                                      const int state_comp, 
                                      const Array4<Real> &forces, 
                                      const int force_comp, const int ncomp,
                                      const Array4<Real> &divu_cc,
                                      const Array4<const Real> &u_mac, 
                                      const Array4<const Real> &v_mac, 
                                      const Array4<const Real> &w_mac,
                                      const Gpu::ManagedVector<BCRec> &BCs,
                                      GpuArray<int,3> const& iconserv,
                                      bool return_state_not_flux)
#endif
