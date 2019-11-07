#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#endif

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MacProjector.H>

#include <incflo.H>
#include <incflo_proj_F.H>

using namespace amrex;

//
// Computes the following decomposition:
// 
//    u + c*grad(phi)/ro = u*  with  div(ep*u) = 0
//
// Inputs:
// 
//   lev               = the AMR level
//   u_mac,v_mac,w_mac = the MAC velocity field to be projected
//   density_in        = the cell-centered density
//
// Outputs:
//
//  phi               = the projection auxiliary function
//  u_mac,v_mac,w_mac = the PROJECTED MAC velocity field 
//
// Notes:
//
//  phi is computed by solving
//
//       div(ep*grad(phi)/rho) = div(ep * u*)
//
//  WARNING: this method returns the MAC velocity with up-to-date BCs in place
// 
void 
incflo::apply_MAC_projection (Vector< std::unique_ptr<MultiFab> >& u_mac, 
                              Vector< std::unique_ptr<MultiFab> >& v_mac,
                              Vector< std::unique_ptr<MultiFab> >& w_mac,
                              Vector< std::unique_ptr<MultiFab> >& density_in,
                              Real time)
{
   BL_PROFILE("incflo::apply_MAC_projection()");

   if (incflo_verbose > 0)
      Print() << "MAC Projection:\n";

   // Check that everything is consistent with amrcore
   // update_internals();

   // Setup for solve
   Vector<Array<MultiFab*,AMREX_SPACEDIM> > mac_vel;
   mac_vel.resize(finest_level+1);

   if (incflo_verbose > 0)
      Print() << " >> Before projection\n" ; 

   // This will hold (1/rho) on faces
   Vector< Array< std::unique_ptr<MultiFab>, AMREX_SPACEDIM> > rho_face;
   rho_face.resize(finest_level+1);

   // This will hold the RHS going into the MAC projection
   Vector< std::unique_ptr<MultiFab> > mac_rhs;
   mac_rhs.resize(finest_level+1);

   // This will hold the solution from the MAC projection
   Vector< std::unique_ptr<MultiFab> > mac_phi;
   mac_phi.resize(finest_level+1);

   for ( int lev=0; lev <= finest_level; ++lev )
   {
      density_in[lev]->FillBoundary(geom[lev].periodicity());

      // We make these with no ghost cells
#ifdef AMREX_USE_EB
      rho_face[lev][0].reset(new MultiFab(u_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]));
      rho_face[lev][1].reset(new MultiFab(v_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]));
      rho_face[lev][2].reset(new MultiFab(w_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]));

      // This doesn't need ghost cells either
      mac_rhs[lev].reset(new MultiFab(grids[lev],dmap[lev],1,0,MFInfo(),*ebfactory[lev]));
      mac_phi[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]));
#else
      rho_face[lev][0].reset(new MultiFab(u_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo()));
      rho_face[lev][1].reset(new MultiFab(v_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo()));
      rho_face[lev][2].reset(new MultiFab(w_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo()));

      // This doesn't need ghost cells either
      mac_rhs[lev].reset(new MultiFab(grids[lev],dmap[lev],1,0,MFInfo()));
      mac_phi[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost,MFInfo()));
#endif

      mac_rhs[lev] -> setVal(0.);
      mac_phi[lev] -> setVal(0.);

      // Define ep and rho on faces
      average_cellcenter_to_face( GetArrOfPtrs(rho_face[lev]), *density_in[lev], geom[lev] );

      int ng_for_invert = 0;
      rho_face[lev][0]->invert(1.0,ng_for_invert);
      rho_face[lev][1]->invert(1.0,ng_for_invert);
      rho_face[lev][2]->invert(1.0,ng_for_invert);

      // Store u_mac in temporaries
      (mac_vel[lev])[0] = u_mac[lev].get();
      (mac_vel[lev])[1] = v_mac[lev].get();
      (mac_vel[lev])[2] = w_mac[lev].get();

      for (int i=0; i<AMREX_SPACEDIM; ++i)
         (mac_vel[lev])[i]->FillBoundary( geom[lev].periodicity() );
      
      // if (incflo_verbose > 0)
      {
#ifdef AMREX_USE_EB
         bool already_on_centroid = true;
         EB_computeDivergence(*mac_rhs[lev],
                              GetArrOfConstPtrs(mac_vel[lev]),
                              geom[lev], already_on_centroid);
#else
         computeDivergence(*mac_rhs[lev],
                           GetArrOfConstPtrs(mac_vel[lev]),
                           geom[lev]);
#endif
         Print() << "  * On level "<< lev
                 << " max(abs(divu)) = " << Norm(mac_rhs,lev,0,0) << "\n";
      }
   }

   Gpu::synchronize();

   //
   // If we want to set max_coarsening_level we have to send it in to the constructor
   //
   LPInfo lp_info;
   lp_info.setMaxCoarseningLevel(mac_mg_max_coarsening_level);

   //
   // Perform MAC projection
   //
   MacProjector macproj( mac_vel, GetVecOfArrOfPtrsConst(rho_face), geom, lp_info);

   // The boundary conditions need only be set at level 0
   int bc_lo[3], bc_hi[3];
   Box domain(geom[0].Domain());
   set_ppe_bcs(bc_lo, bc_hi,
               domain.loVect(), domain.hiVect(),
               &nghost,
               bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
               bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
               bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

   macproj.setDomainBC({(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]},
                       {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]});

   macproj.setVerbose   ( mac_mg_verbose);
   macproj.setCGVerbose ( mac_mg_cg_verbose);
   macproj.setMaxIter   ( mac_mg_maxiter);
   macproj.setCGMaxIter ( mac_mg_cg_maxiter);   
   // The default bottom solver is BiCG
   // Other options include:
   ///   Hypre IJ AMG solver
   //    macproj.getMLMG().setBottomSolver(MLMG::BottomSolver::hypre);
   ///   regular smoothing
   //    macproj.getMLMG().setBottomSolver(MLMG::BottomSolver::smoother);

   if (mac_bottom_solver_type == "smoother")
   {
      macproj.setBottomSolver(MLMG::BottomSolver::smoother);
   }
   else if (mac_bottom_solver_type == "cg")
   {
      macproj.setBottomSolver(MLMG::BottomSolver::cg);
   }
   else if (mac_bottom_solver_type == "bicgcg")
   {
      macproj.setBottomSolver(MLMG::BottomSolver::bicgcg);
   }
   else if (mac_bottom_solver_type == "cgbicg")
   {
      macproj.setBottomSolver(MLMG::BottomSolver::cgbicg);
   }
   else if (mac_bottom_solver_type == "hypre")
   {
      macproj.setBottomSolver(MLMG::BottomSolver::hypre);
   }

   if (steady_state)
   {
       // Solve using mac_phi as an initial guess -- note that mac_phi is
       //       stored from iteration to iteration
       macproj.project(GetVecOfPtrs(mac_phi), mac_mg_rtol,mac_mg_atol,MLMG::Location::FaceCentroid);
   } 
   else 
   {
       // Solve with initial guess of zero
       macproj.project(mac_mg_rtol,mac_mg_atol,MLMG::Location::FaceCentroid);
   }

   if (incflo_verbose > 0)
      Print() << " >> After projection\n" ; 

   for ( int lev=0; lev <= finest_level ; ++lev )
   {   
      if (incflo_verbose > 0)
      {
         mac_vel[lev][0]->FillBoundary( geom[lev].periodicity() );
         mac_vel[lev][1]->FillBoundary( geom[lev].periodicity() );
         mac_vel[lev][2]->FillBoundary( geom[lev].periodicity() );
         
#ifdef AMREX_USE_EB
         bool already_on_centroid = true;
         EB_computeDivergence(*mac_rhs[lev],
                              GetArrOfConstPtrs(mac_vel[lev]),
                              geom[lev], already_on_centroid);
#else
         computeDivergence(*mac_rhs[lev],
                           GetArrOfConstPtrs(mac_vel[lev]),
                           geom[lev]);
#endif

         Print() << "  * On level "<< lev
                 << " max(abs(divu)) = " << Norm(mac_rhs,lev,0,0) << "\n";
      } 

      // Set bcs on u_mac
      set_MAC_velocity_bcs( lev, mac_rhs, u_mac, v_mac, w_mac, time );

      u_mac[lev]->FillBoundary(geom[lev].periodicity());
      v_mac[lev]->FillBoundary(geom[lev].periodicity());
      w_mac[lev]->FillBoundary(geom[lev].periodicity());
   }

   Gpu::synchronize();
}
