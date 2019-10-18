#include <diffusion_F.H>
#include <incflo.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>
#include <AMReX_EBMultiFabUtil.H>

// For multigrid
#include <AMReX_MLTensorOp.H>
#include <AMReX_MLEBTensorOp.H>

//
// Computation of divtau
//
void
incflo::ComputeDivTau (       Vector< std::unique_ptr<MultiFab> >& divtau_out,
                              Vector< std::unique_ptr<MultiFab> >& vel_in,
                        const Vector< std::unique_ptr<MultiFab> >& ro_in,
                        const Vector< std::unique_ptr<MultiFab> >& eta_in)
{
   BL_PROFILE("incflo::ComputeDivTau");

   int nlev = divtau_out.size();

   Vector<std::unique_ptr<MultiFab> > divtau_aux(nlev);
   for (int lev = 0; lev < nlev; lev++)
   {
     divtau_aux[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost,
           MFInfo(), *ebfactory[lev]));
     divtau_aux[lev]->setVal(0.0);
   }

   // Swap ghost cells and apply BCs to velocity
   Real time = 0;
   int extrap_dir_bcs = 0;
   incflo_set_velocity_bcs (time, vel, extrap_dir_bcs);

   // The boundary conditions need only be set once -- we do this at level 0
   int bc_lo[3], bc_hi[3];

   // Whole domain
   Box domain(geom[0].Domain());

   // Set BCs for Poisson's solver
   set_diff_bc (bc_lo, bc_hi,
                domain.loVect(), domain.hiVect(),
                &nghost,
                bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
                bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
                bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

   //
   // First define the operator "ebtensorop"
   //
   //       (alpha * a - beta * (del dot b grad)) sol
   //
   // LPInfo                       info;
   MLEBTensorOp ebtensorop(geom, grids, dmap, LPInfo().setMaxCoarseningLevel(0),
                           amrex::GetVecOfConstPtrs(ebfactory));

   // It is essential that we set MaxOrder of the solver to 2
   // if we want to use the standard sol(i)-sol(i-1) approximation
   // for the gradient at Dirichlet boundaries.
   // The solver's default order is 3 and this uses three points for the
   // gradient at a Dirichlet boundary.
   ebtensorop.setMaxOrder(2);

   // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
   ebtensorop.setDomainBC ( {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]},
                            {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]} );

   // Return div (mu grad)) phi
   ebtensorop.setScalars(0.0, -1.0);

   amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>, AMREX_SPACEDIM>> b;
   b.resize(max_level + 1);

   // Compute the coefficients
   for (int lev = 0; lev < nlev; lev++)
   {
       for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
       {
           BoxArray edge_ba = grids[lev];
           edge_ba.surroundingNodes(dir);
           b[lev][dir].reset(new MultiFab(edge_ba, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));
       }

       average_cellcenter_to_face( GetArrOfPtrs(b[lev]), *eta_in[lev], geom[lev] );

       b[lev][0] -> FillBoundary(geom[lev].periodicity());
       b[lev][1] -> FillBoundary(geom[lev].periodicity());
       b[lev][2] -> FillBoundary(geom[lev].periodicity());

       ebtensorop.setShearViscosity  (lev, GetArrOfConstPtrs(b[lev]));
       ebtensorop.setEBShearViscosity(lev, (*eta[lev]));

       ebtensorop.setLevelBC ( lev, GetVecOfConstPtrs(vel)[lev] );
   }

   MLMG solver(ebtensorop);

   solver.apply(GetVecOfPtrs(divtau_aux), GetVecOfPtrs(vel));

   for (int lev = 0; lev < nlev; lev++)
   {
       incflo_redistribute( lev, *divtau_aux[lev], *divtau_out[lev], 0, AMREX_SPACEDIM );

      // Divide by density
      for (int n = 0; n < 3; n++)
          MultiFab::Divide( *divtau_out[lev], *ro_in[lev], 0, n, 1, 0 );
   }
}
