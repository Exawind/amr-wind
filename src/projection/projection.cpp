#include <AMReX_ParmParse.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_VisMF.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLNodeLaplacian.H>

#include <incflo.H>
#include <mac_F.H>
#include <projection_F.H>

//
// TODO:
// explain that scaling_factor = dt expect when called during initial_projection, when it is set to unity
//
//
// Computes the following decomposition:
//
//    u + grad(phi)/ro = u*,     where div(u) = 0
//
// where u* is a non-div-free velocity field, stored
// by components in u, v, and w. The resulting div-free
// velocity field, u, overwrites the value of u* in u, v, and w.
//
// phi is an auxiliary function related to the pressure p by the relation:
//
//     new p  = phi
//
//     except in the initial iterations when
//
//     new p  = old p + phi
void incflo::ApplyProjection(Real time, Real scaling_factor, bool proj_2)
{
	BL_PROFILE("incflo::ApplyProjection");

    // The boundary conditions need only be set once -- we do this at level 0
	int bc_lo[3], bc_hi[3];

	// Whole domain
    Box domain(geom[0].Domain());

	// Set BCs for Poisson solver
    set_ppe_bc(bc_lo, bc_hi,
               domain.loVect(),
               domain.hiVect(),
               &nghost,
               bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
               bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
               bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

    for(int lev = 0; lev <= finest_level; lev++)
    {
        vel[lev]->FillBoundary(geom[lev].periodicity());
    }

	// Swap ghost cells and apply BCs to velocity
	FillVelocityBC(time, 0);

    ComputeDivU(time);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // Print info about predictor step
        if(incflo_verbose > 1)
        {
            amrex::Print() << "Before projection:" << std::endl;
            PrintMaxValues(time);
        }

        // Here we add the (1/rho gradp) back to ustar (note the +dt)
        if(proj_2)
        {
            // Convert velocities to momenta
            for(int dir = 0; dir < 3; dir++)
            {
                MultiFab::Multiply(*vel[lev], (*ro[lev]), 0, dir, 1, vel[lev]->nGrow());
            }

            MultiFab::Saxpy(*vel[lev], scaling_factor, *gp[lev], 0, 0, 3, vel[lev]->nGrow());

            // Convert momenta back to velocities
            for(int dir = 0; dir < 3; dir++)
            {
                MultiFab::Divide(*vel[lev], (*ro[lev]), 0, dir, 1, vel[lev]->nGrow());
            }
        }
    }

    FillVelocityBC(time, 0);

    // Compute right hand side, AKA div(u)/dt
    ComputeDivU(time);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        divu[lev]->mult(1.0 / scaling_factor, divu[lev]->nGrow());

        // Initialize phi to zero (any non-zero bc's are stored in p0)
        phi[lev]->setVal(0.);
    }

    // Compute the PPE coefficients ( = 1 / rho )
    ComputePoissonCoeff();

	Vector<std::unique_ptr<MultiFab>> fluxes;
    fluxes.resize(max_level + 1);
    for(int lev = 0; lev <= finest_level; lev++)
    {
        fluxes[lev].reset(new MultiFab(vel[lev]->boxArray(),
                                       vel[lev]->DistributionMap(),
                                       vel[lev]->nComp(), 1,
                                       MFInfo(), *ebfactory[lev]));
        fluxes[lev]->setVal(1.0e200);
    }

	// Solve PPE
	SolvePoissonEquation(bcoeff, phi, divu, fluxes, bc_lo, bc_hi);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // NOTE: THE SIGN OF DT (scaling_factor) IS CORRECT HERE
        if(incflo_verbose > 1)
            amrex::Print() << "Multiplying fluxes at level " << lev
                           << " by dt " << scaling_factor << std::endl;

        // The fluxes currently hold MINUS (1/rho) * grad(phi) so we multiply by dt
        fluxes[lev]->mult(scaling_factor, fluxes[lev]->nGrow());

        // Now we correct the velocity with MINUS dt * (1/rho) * grad(phi),
        MultiFab::Add(*vel[lev], *fluxes[lev], 0, 0, 3, 0);

        // The fluxes currently hold MINUS dt * (1/rho) * grad(phi),
        // so now we multiply by rho and divide by (-dt) to get grad(phi)
        fluxes[lev]->mult( -1.0 / scaling_factor, fluxes[lev]->nGrow());
        for(int dir = 0; dir < 3; dir++)
            MultiFab::Multiply(*fluxes[lev], (*ro[lev]), 0, dir, 1, fluxes[lev]->nGrow());

        if(proj_2)
        {
            // p := phi
            MultiFab::Copy(*p[lev], *phi[lev], 0, 0, 1, phi[lev]->nGrow());
            MultiFab::Copy(*gp[lev], *fluxes[lev], 0, 0, 3, fluxes[lev]->nGrow());
        }
        else
        {
            // p := p + phi
            MultiFab::Add(*p[lev], *phi[lev], 0, 0, 1, phi[lev]->nGrow());
            MultiFab::Add(*gp[lev], *fluxes[lev], 0, 0, 3, fluxes[lev]->nGrow());
        }
    }

    // Swap ghost cells and apply BCs to velocity
    FillVelocityBC(time, 0);

    ComputeDivU(time);

    // Print info about predictor step
    if(incflo_verbose > 1)
    {
        amrex::Print() << "After projection: " << std::endl; 
        PrintMaxValues(time);
    }
}

// TODO: See if we can avoid setting everything up in every time step
//       Can we store previous solution as initial guess for next one?
//
// Solve Poisson Equation:
//
//                  div( 1/rho * grad(phi) ) = div(u)
//
void incflo::SolvePoissonEquation(Vector< Vector< std::unique_ptr<MultiFab> > >& rho_inv,
			                      Vector< std::unique_ptr<MultiFab> >& this_phi,
			                      Vector< std::unique_ptr<MultiFab> >& rhs,
			                      Vector< std::unique_ptr<MultiFab> >& fluxes,
			                      int bc_lo[], int bc_hi[])
{
    BL_PROFILE("incflo::SolvePoissonEquation");

    if (nodal_pressure == 1)
    {
        // First define the matrix.
        // Class MLNodeLaplacian describes the following operator:
        //
        //       del dot (sigma grad) phi = rhs,
        //
        // where phi and rhs are nodal, and sigma is cell-centered
        LPInfo info;
        MLNodeLaplacian matrix(geom, grids, dmap, info, amrex::GetVecOfConstPtrs(ebfactory));

        matrix.setGaussSeidel(true);
        matrix.setHarmonicAverage(false);
        matrix.setCoarseningStrategy(MLNodeLaplacian::CoarseningStrategy::Sigma);

        matrix.setDomainBC({(LinOpBCType) bc_lo[0], (LinOpBCType) bc_lo[1], (LinOpBCType) bc_lo[2]},
                           {(LinOpBCType) bc_hi[0], (LinOpBCType) bc_hi[1], (LinOpBCType) bc_hi[2]});

        for (int lev = 0; lev <= finest_level; lev++)
        {
            matrix.setSigma(lev, *(rho_inv[lev][0]));

            // By this point we must have filled the Dirichlet values of phi in ghost cells
            this_phi[lev]->setVal(0.);
            matrix.setLevelBC(lev, GetVecOfConstPtrs(this_phi)[lev]);
        }

        // Then setup the solver ----------------------
        MLMG  solver(matrix);

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

        // This ensures that ghost cells of phi are correctly filled when returned from solver
        solver.setFinalFillBC(true);

        // Finally, solve the system
        solver.solve(GetVecOfPtrs(this_phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol);

        // Get fluxes
        solver.getFluxes(amrex::GetVecOfPtrs(fluxes));
    }
    else
    {
        // First define the matrix.
        // Class MLABecLaplacian describes the following operator:
        //
        //       (alpha * a - beta * (del dot b grad)) phi = rhs
        //
        LPInfo info;
        MLEBABecLap matrix(geom, grids, dmap, info, amrex::GetVecOfConstPtrs(ebfactory));

        // Set alpha and beta
        matrix.setScalars(0.0, -1.0);

        // It is essential that we set MaxOrder of the solver to 2 if we want to use the standard
        // phi(i)-phi(i-1) approximation for the gradient at Dirichlet boundaries.
        // The solver's default order is 3 and this uses three points for the gradient.
        matrix.setMaxOrder(2);

        matrix.setDomainBC({(LinOpBCType) bc_lo[0], (LinOpBCType) bc_lo[1], (LinOpBCType) bc_lo[2]},
                           {(LinOpBCType) bc_hi[0], (LinOpBCType) bc_hi[1], (LinOpBCType) bc_hi[2]});

        for (int lev = 0; lev <= finest_level; lev++)
        {
            // Copy the coefficient into the proper data strutcure
            Vector<const MultiFab*> tmp;
            std::array<MultiFab const*,AMREX_SPACEDIM> b_tmp;

            tmp = amrex::GetVecOfConstPtrs( rho_inv[lev] ) ;
            b_tmp[0] = tmp[0];
            b_tmp[1] = tmp[1];
            b_tmp[2] = tmp[2];

            matrix.setBCoeffs( lev, b_tmp );
            
            // By this point we must have filled the Dirichlet values of phi in ghost cells
            this_phi[lev]->setVal(0.);
            matrix.setLevelBC(lev, GetVecOfConstPtrs(this_phi)[lev]);
        }

        // Then setup the solver ----------------------
        MLMG  solver(matrix);

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

        // This ensures that ghost cells of phi are correctly filled when returned from solver
        solver.setFinalFillBC(true);

        // Finally, solve the system
        solver.solve(GetVecOfPtrs(this_phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol);

        // Get fluxes
        solver.getFluxes(amrex::GetVecOfPtrs(fluxes), MLMG::Location::CellCenter);
    }
    for (int lev = 0; lev <= finest_level; lev++)
    {
       this_phi[lev]->FillBoundary(geom[lev].periodicity());
    }
}


//
// Computes bcoeff = 1/ro at the faces of the scalar cells
//
void incflo::ComputePoissonCoeff()
{
	BL_PROFILE("incflo::ComputePoissonCoeff");

	// Directions
	int xdir = 1;
	int ydir = 2;
	int zdir = 3;

    for(int lev = 0; lev <= finest_level; lev++)
    {
        if (nodal_pressure == 1)
        {
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(*ro[lev],true); mfi.isValid(); ++mfi)
            {
                // Cell-centered tilebox
                Box bx = mfi.tilebox();

                // X direction
                compute_bcoeff_nd(BL_TO_FORTRAN_BOX(bx),
                                  BL_TO_FORTRAN_ANYD((*(bcoeff[lev][0]))[mfi]),
                                  BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
                                  &xdir );

                // Y direction
                compute_bcoeff_nd(BL_TO_FORTRAN_BOX(bx),
                                  BL_TO_FORTRAN_ANYD((*(bcoeff[lev][1]))[mfi]),
                                  BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
                                  &ydir );

                // Z direction
                compute_bcoeff_nd(BL_TO_FORTRAN_BOX(bx),
                                   BL_TO_FORTRAN_ANYD((*(bcoeff[lev][2]))[mfi]),
                                   BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
                                   &zdir );
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel
#endif
            for(MFIter mfi(*ro[lev], true); mfi.isValid(); ++mfi)
            {
                // Tileboxes for staggered components
                Box ubx = mfi.tilebox(e_x);
                Box vbx = mfi.tilebox(e_y);
                Box wbx = mfi.tilebox(e_z);

                // X direction
                compute_bcoeff_mac(BL_TO_FORTRAN_BOX(ubx),
                                   BL_TO_FORTRAN_ANYD((*(bcoeff[lev][0]))[mfi]),
                                   BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
                                   &xdir);

                // Y direction
                compute_bcoeff_mac(BL_TO_FORTRAN_BOX(vbx),
                                   BL_TO_FORTRAN_ANYD((*(bcoeff[lev][1]))[mfi]),
                                   BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
                                   &ydir);

                // Z direction
                compute_bcoeff_mac(BL_TO_FORTRAN_BOX(wbx),
                                   BL_TO_FORTRAN_ANYD((*(bcoeff[lev][2]))[mfi]),
                                   BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
                                   &zdir);
            }
        }

        bcoeff[lev][0]->FillBoundary(geom[lev].periodicity());
        bcoeff[lev][1]->FillBoundary(geom[lev].periodicity());
        bcoeff[lev][2]->FillBoundary(geom[lev].periodicity());
    }
}
