
#include <incflo.H>
#include <bc_mod_F.H>

//
// Set the BCs for face-centroid-based velocity components only
// 
void
incflo::set_MAC_velocity_bcs ( int lev,
                               Vector< std::unique_ptr<MultiFab> >& u_mac,
                               Vector< std::unique_ptr<MultiFab> >& v_mac,
                               Vector< std::unique_ptr<MultiFab> >& w_mac,
                               amrex::Real time)
{
  BL_PROFILE("MacProjection::set_MAC_velocity_bcs()");

  u_mac[lev] -> FillBoundary( geom[lev].periodicity() );
  v_mac[lev] -> FillBoundary( geom[lev].periodicity() );
  w_mac[lev] -> FillBoundary( geom[lev].periodicity() );
    
  Box domain(geom[lev].Domain()); 

  auto lprobtype = probtype;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi((*vel[lev])); mfi.isValid(); ++mfi) // Note we only use vel to define the right MFIter
  {
    const Box& ubx = (*u_mac[lev])[mfi].box();
    IntVect ubx_lo(ubx.loVect());
    IntVect ubx_hi(ubx.hiVect());

    const Box& vbx = (*v_mac[lev])[mfi].box();
    IntVect vbx_lo(vbx.loVect());
    IntVect vbx_hi(vbx.hiVect());

    const Box& wbx = (*w_mac[lev])[mfi].box();
    IntVect wbx_lo(wbx.loVect());
    IntVect wbx_hi(wbx.hiVect());

    IntVect dom_lo(domain.loVect());
    IntVect dom_hi(domain.hiVect());

    Array4<Real> const& u = u_mac[lev]->array(mfi);
    Array4<Real> const& v = v_mac[lev]->array(mfi);
    Array4<Real> const& w = w_mac[lev]->array(mfi);

    IntVect u_lo((*u_mac[lev])[mfi].loVect());
    IntVect u_hi((*u_mac[lev])[mfi].hiVect());

    IntVect v_lo((*v_mac[lev])[mfi].loVect());
    IntVect v_hi((*w_mac[lev])[mfi].hiVect());

    IntVect w_lo((*v_mac[lev])[mfi].loVect());
    IntVect w_hi((*w_mac[lev])[mfi].hiVect());

    Array4<int> const& bct_ilo = bc_ilo[lev]->array();
    Array4<int> const& bct_ihi = bc_ihi[lev]->array();
    Array4<int> const& bct_jlo = bc_jlo[lev]->array();
    Array4<int> const& bct_jhi = bc_jhi[lev]->array();
    Array4<int> const& bct_klo = bc_klo[lev]->array();
    Array4<int> const& bct_khi = bc_khi[lev]->array();

    const int nlft = std::max(0, dom_lo[0]-ubx_lo[0]);
    const int nbot = std::max(0, dom_lo[1]-vbx_lo[1]);
    const int ndwn = std::max(0, dom_lo[2]-wbx_lo[2]);

    const int nrgt = std::max(0, ubx_hi[0]-dom_hi[0]);
    const int ntop = std::max(0, vbx_hi[1]-dom_hi[1]);
    const int nup  = std::max(0, wbx_hi[2]-dom_hi[2]);

    // Create InVects for following Boxes
    IntVect ulo_bx_yz_lo(ubx_lo);
    IntVect ulo_bx_yz_hi(ubx_hi);
    IntVect uhi_bx_yz_lo(ubx_lo);
    IntVect uhi_bx_yz_hi(ubx_hi);

    IntVect vlo_bx_xz_lo(vbx_lo);
    IntVect vlo_bx_xz_hi(vbx_hi);
    IntVect vhi_bx_xz_lo(vbx_lo);
    IntVect vhi_bx_xz_hi(vbx_hi);

    IntVect wlo_bx_xy_lo(wbx_lo);
    IntVect wlo_bx_xy_hi(wbx_hi);
    IntVect whi_bx_xy_lo(wbx_lo);
    IntVect whi_bx_xy_hi(wbx_hi);

    // Fix lo and hi limits
    // Box 'yz'
    ulo_bx_yz_lo[0] = dom_lo[0];
    ulo_bx_yz_hi[0] = dom_lo[0];

    uhi_bx_yz_lo[0] = dom_hi[0]+1;
    uhi_bx_yz_hi[0] = dom_hi[0]+1;

    // Box 'xz'
    vlo_bx_xz_lo[1] = dom_lo[1];
    vlo_bx_xz_hi[1] = dom_lo[1];

    vhi_bx_xz_lo[1] = dom_hi[1]+1;
    vhi_bx_xz_hi[1] = dom_hi[1]+1;

    // Box 'xy'
    wlo_bx_xy_lo[2] = dom_lo[2];
    wlo_bx_xy_hi[2] = dom_lo[2];

    whi_bx_xy_lo[2] = dom_hi[2]+1;
    whi_bx_xy_hi[2] = dom_hi[2]+1;

    // Create 2D boxes for CUDA loops
    const Box ulo_bx_yz(ulo_bx_yz_lo, ulo_bx_yz_hi);
    const Box uhi_bx_yz(uhi_bx_yz_lo, uhi_bx_yz_hi);

    const Box vlo_bx_xz(vlo_bx_xz_lo, vlo_bx_xz_hi);
    const Box vhi_bx_xz(vhi_bx_xz_lo, vhi_bx_xz_hi);

    const Box wlo_bx_xy(wlo_bx_xy_lo, wlo_bx_xy_hi);
    const Box whi_bx_xy(whi_bx_xy_lo, whi_bx_xy_hi);

    const int minf = bc_list.get_minf();
    const int  nsw = bc_list.get_nsw();
    const int slip = bc_list.get_slip();
    const int wall_model = bc_list.get_wall_model();
      
    amrex::Real* p_bc_u = m_bc_u.data();
    amrex::Real* p_bc_v = m_bc_v.data();
    amrex::Real* p_bc_w = m_bc_w.data();

    // NOTE - we only call this for MAC velocities which are only defined on normal faces

    if (nlft > 0)
    {
      amrex::ParallelFor(ulo_bx_yz,
        [bct_ilo,dom_lo,dom_hi,minf,nsw,slip,wall_model,lprobtype,p_bc_u,u] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
        const int bct = bct_ilo(dom_lo[0]-1,j,k,0);
        if (bct == minf) u(i,j,k) = p_bc_u[bcv];
        if (bct ==  nsw || bct == slip || bct == wall_model) u(i,j,k) = 0.;

        if (bct == minf && lprobtype == 31) 
        {
           Real y = (j+0.5) / (dom_hi[1] - dom_lo[1] + 1);
           u(i,j,k) = 6.0 * p_bc_u[bcv] * y * (1.0 - y);
        }
      });
    }

    if (nrgt > 0)
    {
      amrex::ParallelFor(uhi_bx_yz,
        [bct_ihi,dom_lo,dom_hi,minf,nsw,slip,wall_model,lprobtype,p_bc_u,u] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
        const int bct = bct_ihi(dom_hi[0]+1,j,k,0);
        if (bct == minf) u(i,j,k) = p_bc_u[bcv];
        if (bct ==  nsw || bct == slip || bct == wall_model) u(i,j,k) = 0.;
      });
    }

    if (nbot > 0)
    {
      amrex::ParallelFor(vlo_bx_xz,
        [bct_jlo,dom_lo,dom_hi,minf,nsw,slip,wall_model,lprobtype,p_bc_v,v] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
        const int bct = bct_jlo(i,dom_lo[1]-1,k,0);
        if (bct == minf) v(i,j,k) = p_bc_v[bcv];
        if (bct ==  nsw || bct == slip || bct == wall_model) v(i,j,k) = 0.;

        if (bct == minf && lprobtype == 32) 
        {
           Real z = (k+0.5) / (dom_hi[2] - dom_lo[2] + 1);
           v(i,j,k) = 6.0 * p_bc_v[bcv] * z * (1.0 - z);
        }
      });
    }

    if (ntop > 0)
    {
      amrex::ParallelFor(vhi_bx_xz,
        [bct_jhi,dom_lo,dom_hi,minf,nsw,slip,wall_model,lprobtype,p_bc_v,v] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
        const int bct = bct_jhi(i,dom_hi[1]+1,k,0);
        if (bct == minf) v(i,j,k) = p_bc_v[bcv];
        if (bct ==  nsw || bct == slip || bct == wall_model) v(i,j,k) = 0.;
      });
    }

    if (ndwn > 0)
    {
      amrex::ParallelFor(wlo_bx_xy,
        [bct_klo,dom_lo,dom_hi,minf,nsw,slip,wall_model,lprobtype,p_bc_w,w] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
        const int bct = bct_klo(i,j,dom_lo[2]-1,0);
        if (bct == minf) w(i,j,k) = p_bc_w[bcv];
        if (bct ==  nsw || bct == slip || bct == wall_model) w(i,j,k) = 0.;

        if (bct == minf && lprobtype == 33) 
        {
           Real x = (i+0.5) / (dom_hi[0] - dom_lo[0] + 1);
           w(i,j,k) = 6.0 * p_bc_w[bcv] * x * (1.0 - x);
        }
      });
    }

    if (nup > 0)
    {
      amrex::ParallelFor(whi_bx_xy,
        [bct_khi,dom_lo,dom_hi,minf,nsw,slip,wall_model,lprobtype,p_bc_w,w] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
        const int bct = bct_khi(i,j,dom_hi[2]+1,0);
        if (bct == minf) w(i,j,k) = p_bc_w[bcv];
        if (bct ==  nsw || bct == slip || bct == wall_model) w(i,j,k) = 0.;
      });
    }

  }
}
