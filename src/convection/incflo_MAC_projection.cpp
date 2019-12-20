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
//   u_mac,v_mac,w_mac = the MAC velocity field to be projected
//   density           = the cell-centered density
//
// Outputs:
//
//  u_mac,v_mac,w_mac = the PROJECTED MAC velocity field 
//
// Notes:
//
//  phi, the projection auxiliary function, is computed by solving
//
//       div(ep*grad(phi)/rho) = div(ep * u*)
// 
void 
incflo::apply_MAC_projection (Vector<MultiFab>& u_mac,
                              Vector<MultiFab>& v_mac,
                              Vector<MultiFab>& w_mac,
                              Vector<MultiFab*> const& density,
                              Real time)
{
    BL_PROFILE("incflo::apply_MAC_projection()");

    if (incflo_verbose > 0) amrex::Print() << "MAC Projection:\n";

    // This will hold (1/rho) on faces
    Vector<Array<MultiFab ,AMREX_SPACEDIM> > rho_face(finest_level+1);
    Vector<Array<MultiFab*,AMREX_SPACEDIM> > mac_vec(finest_level+1);
    for (int lev=0; lev <= finest_level; ++lev)
    {
        rho_face[lev][0].define(u_mac[lev].boxArray(),dmap[lev],1,0,MFInfo(),Factory(lev));
        rho_face[lev][1].define(v_mac[lev].boxArray(),dmap[lev],1,0,MFInfo(),Factory(lev));
        rho_face[lev][2].define(w_mac[lev].boxArray(),dmap[lev],1,0,MFInfo(),Factory(lev));

        amrex::average_cellcenter_to_face(GetArrOfPtrs(rho_face[lev]), *density[lev], geom[lev]);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            rho_face[lev][idim].invert(1.0, 0);
        }

        mac_vec[lev][0] = &u_mac[lev];
        mac_vec[lev][1] = &v_mac[lev];
        mac_vec[lev][2] = &w_mac[lev];
    }

    //
    // If we want to set max_coarsening_level we have to send it in to the constructor
    //
    LPInfo lp_info;
    lp_info.setMaxCoarseningLevel(mac_mg_max_coarsening_level);

    //
    // Perform MAC projection
    //
    MacProjector macproj(mac_vec, GetVecOfArrOfConstPtrs(rho_face), Geom(0,finest_level), lp_info);

    macproj.setDomainBC(get_projection_bclo(), get_projection_bchi());

    macproj.project(mac_mg_rtol,mac_mg_atol,MLMG::Location::FaceCentroid);

    VisMF::Write(u_mac[0], "umac");
    VisMF::Write(v_mac[0], "vmac");
    VisMF::Write(w_mac[0], "wmac");
    amrex::Abort("xxxxx end of apply_MAC_projection");
}

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
      Print() << " >> Before MAC projection\n";

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

#ifdef AMREX_USE_EB
    const FabFactory<FArrayBox>& factory =  *ebfactory[lev];
#else
    const FabFactory<FArrayBox>& factory = FArrayBoxFactory();
#endif

      density_in[lev]->FillBoundary(geom[lev].periodicity());

      // We make these with no ghost cells
      rho_face[lev][0].reset(new MultiFab(u_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),factory));
      rho_face[lev][1].reset(new MultiFab(v_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),factory));
      rho_face[lev][2].reset(new MultiFab(w_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),factory));

      // mac_rhs doesn't need ghost cells either
      mac_rhs[lev].reset(new MultiFab(grids[lev],dmap[lev],1,0,MFInfo(),factory));

      // mac_phi only needs one ghost cell
      mac_phi[lev].reset(new MultiFab(grids[lev],dmap[lev],1,1,MFInfo(),factory));

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

   //
   // If we want to set max_coarsening_level we have to send it in to the constructor
   //
   LPInfo lp_info;
   lp_info.setMaxCoarseningLevel(mac_mg_max_coarsening_level);

   //
   // Perform MAC projection
   //
   MacProjector macproj( mac_vel, GetVecOfArrOfPtrsConst(rho_face), geom, lp_info);

   macproj.setDomainBC(ppe_lobc, ppe_hibc);

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
      Print() << " >> After  MAC projection\n" ; 

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
      set_MAC_velocity_bcs( lev, u_mac, v_mac, w_mac, time );

      u_mac[lev]->FillBoundary(geom[lev].periodicity());
      v_mac[lev]->FillBoundary(geom[lev].periodicity());
      w_mac[lev]->FillBoundary(geom[lev].periodicity());
   }
}
