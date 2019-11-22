#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

#include <incflo.H>
#include <bc_mod_F.H>

//
//  These subroutines set the BCs for the vel_arr components only.
//

void
incflo::incflo_set_velocity_bcs (Real time,
                                 Vector< std::unique_ptr<MultiFab> > & vel_in,
                                 int extrap_dir_bcs) const
{
  BL_PROFILE("incflo::incflo_set_velocity_bcs()");

  int nlev = vel_in.size();

  for (int lev = 0; lev < nlev; lev++)
  {
     // Set all values outside the domain to covered_val just to avoid use of undefined
     vel_in[lev]->setDomainBndry(covered_val,geom[lev]);

     vel_in[lev] -> FillBoundary (geom[lev].periodicity());
     Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(*vel_in[lev]); mfi.isValid(); ++mfi) {
         set_velocity_bcs(time, lev, (*vel_in[lev])[mfi], domain, extrap_dir_bcs);
     }

#ifdef AMREX_USE_EB
     EB_set_covered(*vel_in[lev], 0, vel_in[lev]->nComp(), vel_in[lev]->nGrow(), covered_val);
#endif

     // Do this after as well as before to pick up terms that got updated in the call above
     vel_in[lev] -> FillBoundary (geom[lev].periodicity());
  }
}

void
incflo::set_velocity_bcs(Real time,
                         const int lev,
                         FArrayBox& vel_fab,
                         const Box& domain,
                         const int extrap_dir_bcs) const
{
  // Must do this because probtype is a member of incflo so can not be "just passed" to the GPU
  auto lprobtype = probtype;

  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<Real> const& vel_arr = vel_fab.array();

  IntVect vel_lo(vel_fab.loVect());
  IntVect vel_hi(vel_fab.hiVect());

  Array4<const int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = bc_klo[lev]->array();
  Array4<const int> const& bct_khi = bc_khi[lev]->array();

  const int nlft = std::max(0, dom_lo[0]-vel_lo[0]);
  const int nbot = std::max(0, dom_lo[1]-vel_lo[1]);
  const int ndwn = std::max(0, dom_lo[2]-vel_lo[2]);

  const int nrgt = std::max(0, vel_hi[0]-dom_hi[0]);
  const int ntop = std::max(0, vel_hi[1]-dom_hi[1]);
  const int nup  = std::max(0, vel_hi[2]-dom_hi[2]);

  // Create InVects for following 2D Boxes
  IntVect bx_yz_lo_lo_2D(vel_lo), bx_yz_lo_hi_2D(vel_hi);
  IntVect bx_yz_hi_lo_2D(vel_lo), bx_yz_hi_hi_2D(vel_hi);
  IntVect bx_xz_lo_lo_2D(vel_lo), bx_xz_lo_hi_2D(vel_hi);
  IntVect bx_xz_hi_lo_2D(vel_lo), bx_xz_hi_hi_2D(vel_hi);
  IntVect bx_xy_lo_lo_2D(vel_lo), bx_xy_lo_hi_2D(vel_hi);
  IntVect bx_xy_hi_lo_2D(vel_lo), bx_xy_hi_hi_2D(vel_hi);

  // Fix lo and hi limits
  bx_yz_lo_lo_2D[0] = dom_lo[0]-1;
  bx_yz_lo_hi_2D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_2D[0] = dom_hi[0]+1;
  bx_yz_hi_hi_2D[0] = dom_hi[0]+1;

  bx_xz_lo_lo_2D[1] = dom_lo[1]-1;
  bx_xz_lo_hi_2D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_2D[1] = dom_hi[1]+1;
  bx_xz_hi_hi_2D[1] = dom_hi[1]+1;

  bx_xy_lo_lo_2D[2] = dom_lo[2]-1;
  bx_xy_lo_hi_2D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_2D[2] = dom_hi[2]+1;
  bx_xy_hi_hi_2D[2] = dom_hi[2]+1;

  // Create 2D boxes for CUDA loops
  const Box bx_yz_lo_2D(bx_yz_lo_lo_2D, bx_yz_lo_hi_2D);
  const Box bx_yz_hi_2D(bx_yz_hi_lo_2D, bx_yz_hi_hi_2D);

  const Box bx_xz_lo_2D(bx_xz_lo_lo_2D, bx_xz_lo_hi_2D);
  const Box bx_xz_hi_2D(bx_xz_hi_lo_2D, bx_xz_hi_hi_2D);

  const Box bx_xy_lo_2D(bx_xy_lo_lo_2D, bx_xy_lo_hi_2D);
  const Box bx_xy_hi_2D(bx_xy_hi_lo_2D, bx_xy_hi_hi_2D);

  // Create InVects for following 3D Boxes
  IntVect bx_yz_lo_hi_3D(vel_hi), bx_xz_lo_hi_3D(vel_hi), bx_xy_lo_hi_3D(vel_hi);
  IntVect bx_yz_hi_lo_3D(vel_lo), bx_xz_hi_lo_3D(vel_lo), bx_xy_hi_lo_3D(vel_lo);

  // Fix lo and hi limits
  bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

  bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

  bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

  // Create 3D boxes for CUDA loops
  const Box bx_yz_lo_3D(vel_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, vel_hi);

  const Box bx_xz_lo_3D(vel_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, vel_hi);

  const Box bx_xy_lo_3D(vel_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, vel_hi);

  for(unsigned i(1); i <= 6; ++i)
  {
    m_bc_u[i] = get_bc_u(i);
    m_bc_v[i] = get_bc_v(i);
    m_bc_w[i] = get_bc_w(i);

    m_bc_r[i] = get_bc_r(i);
    m_bc_t[i] = get_bc_t(i);
  }

  const int minf = bc_list.get_minf();
  const int pinf = bc_list.get_pinf();
  const int pout = bc_list.get_pout();
  const int  nsw = bc_list.get_nsw();
  const int slip = bc_list.get_slip();
  const int wall_model = bc_list.get_wall_model();

  const amrex::Real* p_bc_u = m_bc_u.data();
  const amrex::Real* p_bc_v = m_bc_v.data();
  const amrex::Real* p_bc_w = m_bc_w.data();

  const Real vx_mean_ground_ = vx_mean_ground;
  const Real vy_mean_ground_ = vy_mean_ground;
  const Real nu_mean_ground_ = nu_mean_ground;
  const Real utau_ = utau;
    
  // Coefficients for linear extrapolation to ghost cells -- divide by 3 below
  Real c0 =  6.0;
  Real c1 = -3.0;
  Real c2 =  0.0;

  // Coefficients for quadratic extrapolation to ghost cells -- divide by 3 below
  // (Comment out to stay linear)
  c0 =  8.0;
  c1 = -6.0;
  c2 =  1.0;

  if (nlft > 0)
  {
    amrex::ParallelFor(bx_yz_lo_3D, 3,
      [bct_ilo,dom_lo,dom_hi,pinf,pout,minf,lprobtype,p_bc_u,vel_arr] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if((bct == pinf) or (bct == pout)) {
        vel_arr(i,j,k,n) = vel_arr(dom_lo[0],j,k,n);

      } else if(bct == minf) {
        if (n == 0)
        {
           vel_arr(i,j,k,n) = p_bc_u[bcv];
           if (lprobtype == 31)
           {
               Real y = (j + 0.5) / (dom_hi[1] - dom_lo[1] + 1);
               vel_arr(i,j,k,n) =  6.0 * p_bc_u[bcv] * y * (1.0 - y);
           }
        }
        else
          vel_arr(i,j,k,n) = 0;
      }
    });
  }

  if (nrgt > 0)
  {
    amrex::ParallelFor(bx_yz_hi_3D, 3,
      [bct_ihi,dom_lo,dom_hi,pinf,pout,minf,lprobtype,p_bc_u,vel_arr] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if((bct == pinf) or (bct == pout)) {
        vel_arr(i,j,k,n) = vel_arr(dom_hi[0],j,k,n);
      } else if(bct == minf)
      {
        if(n == 0)
          vel_arr(i,j,k,n) = p_bc_u[bcv];
        else
          vel_arr(i,j,k,n) = 0;
      }
    });
  }

  if (nbot > 0)
  {
    amrex::ParallelFor(bx_xz_lo_3D, 3,
      [bct_jlo,dom_lo,dom_hi,pinf,pout,minf,lprobtype,p_bc_v,vel_arr] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if((bct == pinf) or (bct == pout)) {
        vel_arr(i,j,k,n) = vel_arr(i,dom_lo[1],k,n);
      } else if(bct == minf)
      {
        if(n == 1)
        {
           vel_arr(i,j,k,n) = p_bc_v[bcv];
           if (lprobtype == 32)
           {
               Real z = (k + 0.5) / (dom_hi[2] - dom_lo[2] + 1);
               vel_arr(i,j,k,n) =  6.0 * p_bc_v[bcv] * z * (1.0 - z);
           }
        }
        else
        {
          vel_arr(i,j,k,n) = 0;
        }
      }
    });
  }

  if (ntop > 0)
  {
    amrex::ParallelFor(bx_xz_hi_3D, 3,
      [bct_jhi,dom_lo,dom_hi,pinf,pout,minf,lprobtype,p_bc_v,vel_arr] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if((bct == pinf) or (bct == pout)) {
        vel_arr(i,j,k,n) = vel_arr(i,dom_hi[1],k,n);
      } else if(bct == minf)
      {
        if(n == 1)
          vel_arr(i,j,k,n) = p_bc_v[bcv];
        else
          vel_arr(i,j,k,n) = 0;
      }
    });
  }

  if (ndwn > 0)
  {
    amrex::ParallelFor(bx_xy_lo_3D, 3,
      [bct_klo,dom_lo,dom_hi,pinf,pout,minf,lprobtype,p_bc_w,vel_arr] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if((bct == pinf) or (bct == pout)) {
        vel_arr(i,j,k,n) = vel_arr(i,j,dom_lo[2],n);
      } else if(bct == minf)
      {
        if(n == 2)
        {
           vel_arr(i,j,k,n) = p_bc_w[bcv];
           if (lprobtype == 33)
           {
               Real x = (i + 0.5) / (dom_hi[0] - dom_lo[0] + 1);
               vel_arr(i,j,k,n) =  6.0 * p_bc_w[bcv] * x * (1.0 - x);
           }
        }
        else
          vel_arr(i,j,k,n) = 0;
      }
    });
  }

  if (nup > 0)
  {
    amrex::ParallelFor(bx_xy_hi_3D, 3,
      [bct_khi,dom_lo,dom_hi,pinf,pout,minf,lprobtype,p_bc_w,vel_arr] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if((bct == pinf) or (bct == pout)) {
        vel_arr(i,j,k,n) = vel_arr(i,j,dom_hi[2],n);
      } else if(bct == minf)
      {
        if(n == 2)
          vel_arr(i,j,k,n) = p_bc_w[bcv];
        else
          vel_arr(i,j,k,n) = 0;
      }
    });
  }

/* *******************************************************************
   Do this section next to make sure nsw over-rides any previous minf
   ****************************************************************** */

  if (nlft > 0)
  {
    amrex::ParallelFor(bx_yz_lo_3D,
      [bct_ilo,dom_lo,nsw,slip,wall_model,p_bc_v,p_bc_w,vel_arr] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if(bct == nsw) {
        vel_arr(i,j,k,0) = 0.;
        vel_arr(i,j,k,1) = p_bc_v[bcv];
        vel_arr(i,j,k,2) = p_bc_w[bcv];
      } else if (bct == slip) {
          vel_arr(i,j,k,0) = 0.0;
          vel_arr(i,j,k,1) = vel_arr(dom_lo[0],j,k,1);
          vel_arr(i,j,k,2) = vel_arr(dom_lo[0],j,k,2);
      } else if (bct == wall_model) {
          amrex::Abort("not implemented for the left wall yet");
      }
    });

    if(extrap_dir_bcs > 0)
    {
      amrex::ParallelFor(bx_yz_lo_2D, 3,
        [bct_ilo,dom_lo,minf,nsw,slip,wall_model,c0,c1,c2,vel_arr] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
        const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

        if(bct == minf || bct == nsw)
          vel_arr(i,j,k,n) = (c0*vel_arr(i,j,k,n) + c1*vel_arr(i+1,j,k,n) + c2*vel_arr(i+2,j,k,n)) / 3.0;

      });
    }
  }

  if (nrgt > 0)
  {
    amrex::ParallelFor(bx_yz_hi_3D,
      [bct_ihi,dom_hi,nsw,slip,wall_model,p_bc_v,p_bc_w,vel_arr] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if(bct == nsw) {
        vel_arr(i,j,k,0) = 0.;
        vel_arr(i,j,k,1) = p_bc_v[bcv];
        vel_arr(i,j,k,2) = p_bc_w[bcv];
      } else if (bct == slip) {
          vel_arr(i,j,k,0) = 0.0;
          vel_arr(i,j,k,1) = vel_arr(dom_hi[0],j,k,1);
          vel_arr(i,j,k,2) = vel_arr(dom_hi[0],j,k,2);
      } else if (bct == wall_model) {
          amrex::Abort("not implemented for the right wall yet");
      }
    });

    if(extrap_dir_bcs > 0)
    {
      amrex::ParallelFor(bx_yz_hi_2D,
        [bct_ihi,dom_hi,minf,nsw,slip,wall_model,c0,c1,c2,vel_arr] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

        if(bct == minf || bct == nsw)
        {
          vel_arr(i,j,k,0) = (c0*vel_arr(i,j,k,0) + c1*vel_arr(i-1,j,k,0) + c2*vel_arr(i-2,j,k,0)) / 3.0;
          vel_arr(i,j,k,1) = (c0*vel_arr(i,j,k,1) + c1*vel_arr(i-1,j,k,1) + c2*vel_arr(i-2,j,k,1)) / 3.0;
          vel_arr(i,j,k,2) = (c0*vel_arr(i,j,k,2) + c1*vel_arr(i-1,j,k,2) + c2*vel_arr(i-2,j,k,2)) / 3.0;
        }
      });
    }
  }

  if (nbot > 0)
  {
    amrex::ParallelFor(bx_xz_lo_3D,
      [bct_jlo,dom_lo,nsw,slip,wall_model,p_bc_u,p_bc_w,vel_arr] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if (bct == nsw) {
        vel_arr(i,j,k,0) = p_bc_u[bcv];
        vel_arr(i,j,k,1) = 0.;
        vel_arr(i,j,k,2) = p_bc_w[bcv];
      } else if (bct == slip) {
          vel_arr(i,j,k,0) = vel_arr(i,dom_lo[1],k,0);
          vel_arr(i,j,k,1) = 0.0;
          vel_arr(i,j,k,2) = vel_arr(i,dom_lo[1],k,2);
      } else if (bct == wall_model) {
          amrex::Abort("not implemented for the bottom wall yet");
      }
    });

    if(extrap_dir_bcs > 0)
    {
      amrex::ParallelFor(bx_xz_lo_2D,
        [bct_jlo,dom_lo,minf,nsw,slip,wall_model,c0,c1,c2,vel_arr] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

        if(bct == minf || bct == nsw)
        {
          vel_arr(i,j,k,0) = (c0*vel_arr(i,j,k,0) + c1*vel_arr(i,j+1,k,0) + c2*vel_arr(i,j+2,k,0)) / 3.0;
          vel_arr(i,j,k,1) = (c0*vel_arr(i,j,k,1) + c1*vel_arr(i,j+1,k,1) + c2*vel_arr(i,j+2,k,1)) / 3.0;
          vel_arr(i,j,k,2) = (c0*vel_arr(i,j,k,2) + c1*vel_arr(i,j+1,k,2) + c2*vel_arr(i,j+2,k,2)) / 3.0;
        }
      });
    }
  }

  if (ntop > 0)
  {
    amrex::ParallelFor(bx_xz_hi_3D,
      [bct_jhi,dom_hi,nsw,slip,wall_model,p_bc_u,p_bc_w,vel_arr] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if (bct == nsw) {
        vel_arr(i,j,k,0) = p_bc_u[bcv];
        vel_arr(i,j,k,1) = 0.;
        vel_arr(i,j,k,2) = p_bc_w[bcv];
      } else if (bct == slip) {
          vel_arr(i,j,k,0) = vel_arr(i,dom_hi[1],k,0);
          vel_arr(i,j,k,1) = 0.0;
          vel_arr(i,j,k,2) = vel_arr(i,dom_hi[1],k,2);
      } else if (bct == wall_model) {
          amrex::Abort("not implemented for the top wall yet");
      }
    });

    if(extrap_dir_bcs > 0)
    {
      amrex::ParallelFor(bx_xz_hi_2D,
        [bct_jhi,dom_hi,minf,nsw,slip,wall_model,c0,c1,c2,vel_arr] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

        if(bct == minf || bct == nsw)
        {
          vel_arr(i,j,k,0) = (c0*vel_arr(i,j,k,0) + c1*vel_arr(i,j-1,k,0) + c2*vel_arr(i,j-2,k,0)) / 3.0 ;
          vel_arr(i,j,k,1) = (c0*vel_arr(i,j,k,1) + c1*vel_arr(i,j-1,k,1) + c2*vel_arr(i,j-2,k,1)) / 3.0 ;
          vel_arr(i,j,k,2) = (c0*vel_arr(i,j,k,2) + c1*vel_arr(i,j-1,k,2) + c2*vel_arr(i,j-2,k,1)) / 3.0 ;
        }
      });
    }
  }

  if (ndwn > 0)
  {
    amrex::ParallelFor(bx_xy_lo_3D,
      [bct_klo,dom_lo,nsw,slip,wall_model,p_bc_u,p_bc_v,vel_arr,vx_mean_ground_,vy_mean_ground_,utau_,nu_mean_ground_] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if (bct == nsw) {
        vel_arr(i,j,k,0) = p_bc_u[bcv];
        vel_arr(i,j,k,1) = p_bc_v[bcv];
        vel_arr(i,j,k,2) = 0.;
      } else if (bct == slip) {
          vel_arr(i,j,k,0) = vel_arr(i,j,dom_lo[2],0);
          vel_arr(i,j,k,1) = vel_arr(i,j,dom_lo[2],1);
          vel_arr(i,j,k,2) = 0.0;
      } else if (bct == wall_model) {
          
          const Real vx = vel_arr(i,j,dom_lo[2],0);
          const Real vy = vel_arr(i,j,dom_lo[2],1);

          // fixme remove this later when we trust this
          const Real vmag = sqrt(pow(vx,2) + pow(vy,2)) + 1.0e-12;
          const Real uh = sqrt(pow(vx_mean_ground_,2) + pow(vy_mean_ground_,2));
          if(vmag/uh > 10.0) printf("uh oh local velocity is large i %d j %d vmag %f uh %f\n",i,j,vmag,uh);
          if(uh/vmag > 10.0) printf("uh oh average velocity is large i %d j %d vmag %f uh %f\n",i,j,vmag,uh);
          
          // simple shear stress model for neutral BL
          // apply as an inhomogeneous Neumann BC
          // fixme this should be the local (nu+nut) but there is a circular dependency on strain rate and bc's
          // we could locally calculate a strain rate but not sure what do with 1 eqn sgs
          vel_arr(i,j,k,0) = utau_*utau_*vx/uh/nu_mean_ground_;
          vel_arr(i,j,k,1) = utau_*utau_*vy/uh/nu_mean_ground_;
          vel_arr(i,j,k,2) = 0.0;
      }
    });

    if(extrap_dir_bcs > 0)
    {
      amrex::ParallelFor(bx_xy_lo_2D,
        [bct_klo,dom_lo,minf,nsw,slip,wall_model,c0,c1,c2,vel_arr] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_klo(i,j,dom_lo[2]-1,0);

        if(bct == minf || bct == nsw)
        {
          vel_arr(i,j,k,0) = (c0*vel_arr(i,j,k,0) + c1*vel_arr(i,j,k+1,0) + c2*vel_arr(i,j,k+2,0)) / 3.0 ;
          vel_arr(i,j,k,1) = (c0*vel_arr(i,j,k,1) + c1*vel_arr(i,j,k+1,1) + c2*vel_arr(i,j,k+2,1)) / 3.0 ;
          vel_arr(i,j,k,2) = (c0*vel_arr(i,j,k,2) + c1*vel_arr(i,j,k+1,2) + c2*vel_arr(i,j,k+2,2)) / 3.0 ;
        }
      });
    }
  }

  if (nup > 0)
  {
    amrex::ParallelFor(bx_xy_hi_3D,
      [bct_khi,dom_hi,nsw,slip,wall_model,p_bc_u,p_bc_v,vel_arr] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if (bct == nsw) {
        vel_arr(i,j,k,0) = p_bc_u[bcv];
        vel_arr(i,j,k,1) = p_bc_v[bcv];
        vel_arr(i,j,k,2) = 0.;
      } else if (bct == slip) {
        vel_arr(i,j,k,0) = vel_arr(i,j,dom_hi[2],0);
        vel_arr(i,j,k,1) = vel_arr(i,j,dom_hi[2],1);
        vel_arr(i,j,k,2) = 0.0;

      } else if (bct == wall_model) {
          amrex::Abort("not implemented for the upper wall yet");
      }
    });

    if(extrap_dir_bcs > 0)
    {
      amrex::ParallelFor(bx_xy_hi_2D,
        [bct_khi,dom_hi,minf,nsw,c0,c1,c2,vel_arr] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_khi(i,j,dom_hi[2]+1,0);

        if(bct == minf || bct == nsw)
        {
          vel_arr(i,j,k,0) = (c0*vel_arr(i,j,k,0) + c1*vel_arr(i,j,k-1,0) + c2*vel_arr(i,j,k-2,0)) / 3.0;
          vel_arr(i,j,k,1) = (c0*vel_arr(i,j,k,1) + c1*vel_arr(i,j,k-1,1) + c2*vel_arr(i,j,k-2,1)) / 3.0;
          vel_arr(i,j,k,2) = (c0*vel_arr(i,j,k,2) + c1*vel_arr(i,j,k-1,2) + c2*vel_arr(i,j,k-2,2)) / 3.0;
        }
      });
    }
  }
}
