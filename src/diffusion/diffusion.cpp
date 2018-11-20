#include <AMReX_ParmParse.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>
#include <diffusion_F.H>
#include <incflo.H>

// For multigrid
#include <AMReX_MLMG.H>
#include <AMReX_MLEBABecLap.H>

//
// Divergence of stress tensor
//
void incflo::ComputeDivTau(int lev,
                                   MultiFab& divtau,
                                   Vector<std::unique_ptr<MultiFab>>& vel_in)
{
	BL_PROFILE("incflo::ComputeDivTau");
	Box domain(geom[lev].Domain());

   // Get EB geometric info
   Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
   Array< const MultiCutFab*,AMREX_SPACEDIM> facecent;
   const amrex::MultiFab*                    volfrac;
   const amrex::MultiCutFab*                 bndrycent;

   areafrac  =   ebfactory[lev] -> getAreaFrac();
   facecent  =   ebfactory[lev] -> getFaceCent();
   volfrac   = &(ebfactory[lev] -> getVolFrac());
   bndrycent = &(ebfactory[lev] -> getBndryCent());

#ifdef _OPENMP
#pragma omp parallel
#endif
   for (MFIter mfi(*vel_in[lev],true); mfi.isValid(); ++mfi) {

      // Tilebox
      Box bx = mfi.tilebox ();

      // this is to check efficiently if this tile contains any eb stuff
      const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel_in[lev])[mfi]);
      const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

      if (flags.getType(bx) == FabType::covered)
      {
         divtau[mfi].setVal(1.2345e200, bx, 0, 3);
      }
      else
      {
         if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular)
         {
            compute_divtau(
               BL_TO_FORTRAN_BOX(bx),
               BL_TO_FORTRAN_ANYD(divtau[mfi]),
               BL_TO_FORTRAN_ANYD((*vel_in[lev])[mfi]),
               (*eta[lev])[mfi].dataPtr(),
               BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
               domain.loVect (), domain.hiVect (),
               bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
               bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
               bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
               geom[lev].CellSize(), &nghost, &explicit_diffusion);
         }
         else
         {
            compute_divtau_eb(
               BL_TO_FORTRAN_BOX(bx),
               BL_TO_FORTRAN_ANYD(divtau[mfi]),
               BL_TO_FORTRAN_ANYD((*vel_in[lev])[mfi]),
               (*eta[lev])[mfi].dataPtr(),
               BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
               BL_TO_FORTRAN_ANYD(flags),
               BL_TO_FORTRAN_ANYD((*areafrac[0])[mfi]),
               BL_TO_FORTRAN_ANYD((*areafrac[1])[mfi]),
               BL_TO_FORTRAN_ANYD((*areafrac[2])[mfi]),
               BL_TO_FORTRAN_ANYD((*facecent[0])[mfi]),
               BL_TO_FORTRAN_ANYD((*facecent[1])[mfi]),
               BL_TO_FORTRAN_ANYD((*facecent[2])[mfi]),
               BL_TO_FORTRAN_ANYD((*volfrac)[mfi]),
               BL_TO_FORTRAN_ANYD((*bndrycent)[mfi]),
               domain.loVect (), domain.hiVect (),
               bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
               bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
               bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
               geom[lev].CellSize(), &nghost, &explicit_diffusion);

         }
      }
   }
}







//
// Implicit diffusion
//
void incflo::DiffuseVelocity(amrex::Real time)
{
	BL_PROFILE("incflo::DiffuseVelocity");

	// Swap ghost cells and apply BCs to velocity
	FillVelocityBC(time, 0);

    // The boundary conditions need only be set once -- we do this at level 0
	int bc_lo[3], bc_hi[3];

	// Whole domain
	Box domain(geom[0].Domain());

	// Set BCs for Poisson solver
    set_diff_bc(bc_lo, bc_hi,
                domain.loVect(), domain.hiVect(),
				&nghost,
                bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
                bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
                bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

    // Declare, resize and reset local variables for matrix solve
    Vector<Vector<std::unique_ptr<MultiFab>>> beta;
    Vector<std::unique_ptr<MultiFab>> sol;
    Vector<std::unique_ptr<MultiFab>> RHS;
    beta.resize(finest_level + 1);
    sol.resize(finest_level + 1);
    RHS.resize(finest_level + 1);
    for(int lev = 0; lev <= finest_level; lev++)
    {
        beta[lev].resize(3);
        for(int dir = 0; dir < 3; dir++)
        {
            BoxArray edge_ba = grids[lev].surroundingNodes(dir);
            beta[lev][dir].reset(new MultiFab(edge_ba, dmap[lev], 1, nghost));
        }
        sol[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
        RHS[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
    }

	// Compute the coefficients
	ComputeDiffusionCoeff(beta);

	// Loop over the velocity components
	for(int dir = 0; dir < 3; dir++)
	{
        if(incflo_verbose > 1)
        {
            amrex::Print() << "Diffusing velocity component " << dir << std::endl;
        }

        for(int lev = 0; lev <= finest_level; lev++)
        {
            RHS[lev]->copy(*vel[lev], dir, 0, 1, nghost, nghost);
            sol[lev]->copy(*vel[lev], dir, 0, 1, nghost, nghost);
        }


		// Solve (1 - div beta grad) u_new = RHS
		// Here RHS = "vel" which is the current approximation to the new-time velocity (without diffusion terms)

		SolveDiffusionEquation(beta, sol, RHS, bc_lo, bc_hi);

        for(int lev = 0; lev <= finest_level; lev++)
        {
            vel[lev]->copy(*sol[lev], 0, dir, 1, nghost, nghost);
        }
	}

	// Swap ghost cells and apply BCs to velocity
	FillVelocityBC(time, 0);
}

// TODO: See if we can avoid setting everything up in every time step
//       Can we store previous solution as initial guess for next one?
//
// Solve :
//                  (1 - div dot eta grad) u = rhs
//
void incflo::SolveDiffusionEquation(Vector<Vector<std::unique_ptr<MultiFab>>>& b,
                                      Vector<std::unique_ptr<MultiFab>>& sol,
                                      Vector<std::unique_ptr<MultiFab>>& RHS,
                                      int bc_lo[], int bc_hi[])
{
	BL_PROFILE("incflo::SolveDiffusionEquation");

    int debug = 0;
	// First define the matrix.
	LPInfo info;
	// Class MLABecLaplacian describes the following operator:
	//
	//       (alpha * a - beta * (del dot b grad)) sol = rhs
	//
    MLEBABecLap matrix(geom, grids, dmap, info, amrex::GetVecOfConstPtrs(ebfactory));

	// Set alpha and beta
	matrix.setScalars(1.0, dt);

    // It is essential that we set MaxOrder of the solver to 2 if we want to use the standard
    // phi(i)-phi(i-1) approximation for the gradient at Dirichlet boundaries.
    // The solver's default order is 3 and this uses three points for the gradient.
	matrix.setMaxOrder(2);

	// LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
	matrix.setDomainBC({(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]},
					   {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]});

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // Copy the PPE coefficient into the proper data strutcure
        Vector<const MultiFab*> tmp;
        std::array<MultiFab const*, AMREX_SPACEDIM> b_tmp;

        tmp = amrex::GetVecOfConstPtrs(b[lev]);
        b_tmp[0] = tmp[0];
        b_tmp[1] = tmp[1];
        b_tmp[2] = tmp[2];

        // Define RHS = (ro) * (vel)
        MultiFab::Multiply((*RHS[lev]), (*ro[lev]), 0, 0, 1, RHS[lev]->nGrow());

        // This sets the spatially varying A coefficients
        matrix.setACoeffs(lev, (*ro[lev]));

        // This sets the spatially varying b coefficients
        matrix.setBCoeffs(lev, b_tmp);

        // By this point we must have filled the Dirichlet values of sol stored in the ghost cells
        amrex::Print() << ++debug << std::endl; 
        matrix.setLevelBC(lev, GetVecOfConstPtrs(sol)[lev]);
        amrex::Print() << ++debug << std::endl; 
    }

	// Then setup the solver ----------------------
	MLMG solver(matrix);
    amrex::Print() << ++debug << std::endl; 

    // The default bottom solver is BiCG
    if(bottom_solver_type == "smoother")
    {
       solver.setBottomSolver(MLMG::BottomSolver::smoother);
    }
    else if(bottom_solver_type == "hypre")
    {
       solver.setBottomSolver(MLMG::BottomSolver::hypre);
    }

	solver.setMaxIter(mg_max_iter);
	solver.setMaxFmgIter(mg_max_fmg_iter);
	solver.setVerbose(mg_verbose);
	solver.setCGVerbose(mg_cg_verbose);
	solver.setCGMaxIter(mg_cg_maxiter);

	// This ensures that ghost cells of sol are correctly filled when returned from the solver
	solver.setFinalFillBC(true);

    amrex::Print() << ++debug << std::endl; 
	// Finally, solve the system
	solver.solve(GetVecOfPtrs(sol), GetVecOfConstPtrs(RHS), mg_rtol, mg_atol);
    amrex::Print() << ++debug << std::endl; 

    for(int lev = 0; lev <= finest_level; lev++)
    {
        sol[lev]->FillBoundary(geom[lev].periodicity());
    }
}

//
// Computes beta = eta at the faces of the scalar cells
//
void incflo::ComputeDiffusionCoeff(Vector<Vector<std::unique_ptr<MultiFab>>>& beta)
{
	BL_PROFILE("incflo::ComputeDiffusionCoeff");

	// Directions
	int xdir = 1;
	int ydir = 2;
	int zdir = 3;

    for(int lev = 0; lev <= finest_level; lev++)
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for(MFIter mfi(*eta[lev], true); mfi.isValid(); ++mfi)
        {
            // Tileboxes for staggered components
            Box ubx = mfi.tilebox(e_x);
            Box vbx = mfi.tilebox(e_y);
            Box wbx = mfi.tilebox(e_z);

            // X direction
            compute_bcoeff_diff(BL_TO_FORTRAN_BOX(ubx),
                                BL_TO_FORTRAN_ANYD((*(beta[lev][0]))[mfi]),
                                BL_TO_FORTRAN_ANYD((*eta[lev])[mfi]),
                                &xdir);

            // Y direction
            compute_bcoeff_diff(BL_TO_FORTRAN_BOX(vbx),
                                BL_TO_FORTRAN_ANYD((*(beta[lev][1]))[mfi]),
                                BL_TO_FORTRAN_ANYD((*eta[lev])[mfi]),
                                &ydir);

            // Z direction
            compute_bcoeff_diff(BL_TO_FORTRAN_BOX(wbx),
                                BL_TO_FORTRAN_ANYD((*(beta[lev][2]))[mfi]),
                                BL_TO_FORTRAN_ANYD((*eta[lev])[mfi]),
                                &zdir);
        }
        beta[lev][0]->FillBoundary(geom[lev].periodicity());
        beta[lev][1]->FillBoundary(geom[lev].periodicity());
        beta[lev][2]->FillBoundary(geom[lev].periodicity());
    }
}
