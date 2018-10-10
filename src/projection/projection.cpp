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
void incflo::incflo_apply_projection(Real time, Real scaling_factor, bool proj_2)
{
	BL_PROFILE("incflo::incflo_apply_projection");

    // Set domain BCs for Poisson solver
    // The domain BCs refer to level 0 only
    int bc_lo[3], bc_hi[3];
    Box domain(geom[0].Domain());

    set_ppe_bc(bc_lo,
               bc_hi,
               domain.loVect(),
               domain.hiVect(),
               &nghost,
               bc_ilo[0]->dataPtr(),
               bc_ihi[0]->dataPtr(),
               bc_jlo[0]->dataPtr(),
               bc_jhi[0]->dataPtr(),
               bc_klo[0]->dataPtr(),
               bc_khi[0]->dataPtr());

    for(int lev = 0; lev < nlev; lev++)
    {
        vel[lev]->FillBoundary(geom[lev].periodicity());
    }

    // Swap ghost cells and apply BCs to velocity
    incflo_set_velocity_bcs(time, 0);

    incflo_compute_divu(time);

    for(int lev = 0; lev < nlev; lev++)
    {
        // Print info about predictor step
        if(verbose > 1)
        {
            amrex::Print() << "At level " << lev << ", before projection: \n";
            incflo_print_max_vel(lev);
            amrex::Print() << "max(abs(divu)) = " << incflo_norm0(divu, lev, 0) << "\n";
        }

        // Here we add the (1/rho gradp) back to ustar (note the +dt)
        if(proj_2)
        {
            // Convert velocities to momenta
            for(int n = 0; n < 3; n++)
                MultiFab::Multiply(*vel[lev], (*ro[lev]), 0, n, 1, vel[lev]->nGrow());

            MultiFab::Saxpy(*vel[lev], scaling_factor, *gp[lev], 0, 0, 3, vel[lev]->nGrow());

            // Convert momenta back to velocities
            for(int n = 0; n < 3; n++)
                MultiFab::Divide(*vel[lev], (*ro[lev]), 0, n, 1, vel[lev]->nGrow());
        }
    }

    incflo_set_velocity_bcs(time, 0);

    // Compute right hand side, AKA div(u)/dt
    incflo_compute_divu(time);

    for(int lev = 0; lev < nlev; lev++)
    {
        divu[lev]->mult(1.0 / scaling_factor, divu[lev]->nGrow());

        // Initialize phi to zero (any non-zero bc's are stored in p0)
        phi[lev]->setVal(0.);
    }

    // Compute the PPE coefficients
    incflo_compute_bcoeff_ppe();

	Vector<std::unique_ptr<MultiFab>> fluxes;
    fluxes.resize(nlev);
    for(int lev = 0; lev < nlev; lev++)
    {
        fluxes[lev].reset(new MultiFab(vel[lev]->boxArray(), 
                                       vel[lev]->DistributionMap(),
                                       vel[lev]->nComp(), 1,
                                       MFInfo(), *ebfactory[lev]));
        fluxes[lev]->setVal(1.0e200);
    }

	// Solve PPE
	solve_poisson_equation(bcoeff, phi, divu, fluxes, bc_lo, bc_hi);

    for(int lev = 0; lev < nlev; lev++)
    {
        // NOTE: THE SIGN OF DT (scaling_factor) IS CORRECT HERE
        if(verbose > 1)
            amrex::Print() << "Multiplying fluxes at level " << lev 
                            << " by dt " << scaling_factor << std::endl;
        fluxes[lev]->mult(scaling_factor, fluxes[lev]->nGrow());
        MultiFab::Add(*vel[lev], *fluxes[lev], 0, 0, 3, 0);

        // After using the fluxes, which currently hold MINUS dt * (1/rho) * grad(phi),
        //    to modify the velocity field,  convert them to hold grad(phi)
        fluxes[lev]->mult(-1.0 / scaling_factor, fluxes[lev]->nGrow());
        for(int n = 0; n < 3; n++)
            MultiFab::Multiply(*fluxes[lev], (*ro[lev]), 0, n, 1, fluxes[lev]->nGrow());

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
    incflo_set_velocity_bcs(time, 0);

    incflo_compute_divu(time);

    // Print info about predictor step
    if(verbose > 1)
    {
        for(int lev = 0; lev < nlev; lev++)
        {
            amrex::Print() << "At level " << lev << ", after projection: \n";
            incflo_print_max_vel(lev);
            amrex::Print() << "max(abs(divu)) = " << incflo_norm0(divu, lev, 0) << "\n";
        }
    }
}

//
// Solve PPE:
//
//                  div( 1/rho * grad(phi) ) = div(u)
//
void incflo::solve_poisson_equation(Vector<Vector<std::unique_ptr<MultiFab>>>& b,
                                    Vector<std::unique_ptr<MultiFab>>& this_phi,
                                    Vector<std::unique_ptr<MultiFab>>& rhs,
                                    Vector<std::unique_ptr<MultiFab>>& fluxes,
                                    int bc_lo[], int bc_hi[])
{
	BL_PROFILE("incflo::solve_poisson_equation");

    //
    // First define the matrix (operator).
    // Class MLABecLaplacian describes the following operator:
    //
    //       (alpha * a - beta * (del dot b grad)) phi
    //
    LPInfo info;
    MLEBABecLap matrix(geom, grids, dmap, info, GetVecOfConstPtrs(ebfactory));
    Vector<const MultiFab*> tmp;
    std::array<MultiFab const*, AMREX_SPACEDIM> b_tmp;

    int lev = 0;

    // Copy the PPE coefficient into the proper data strutcure
    tmp = GetVecOfConstPtrs(b[lev]);
    b_tmp[0] = tmp[0];
    b_tmp[1] = tmp[1];
    b_tmp[2] = tmp[2];

    // It is essential that we set MaxOrder of the solver to 2
    // if we want to use the standard phi(i)-phi(i-1) approximation
    // for the gradient at Dirichlet boundaries.
    // The solver's default order is 3 and this uses three points for the
    // gradient at a Dirichlet boundary.
    matrix.setMaxOrder(2);

    // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
    matrix.setDomainBC({(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]},
                       {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]});

    matrix.setScalars(0.0, -1.0);
    matrix.setBCoeffs(lev, b_tmp);

    // By this point we must have filled the Dirichlet values of phi stored in the ghost cells
    this_phi[lev]->setVal(0.);
    matrix.setLevelBC(lev, GetVecOfConstPtrs(this_phi)[lev]);

    //
    // Then setup the solver ----------------------
    //
    MLMG solver(matrix);

   // The default bottom solver is BiCG
   // Other options include: 
   ///   regular smoothing ("smoother")
   ///   Hypre IJ AMG solver ("hypre")
   if (bottom_solver_type == "smoother")
   { 
      solver.setBottomSolver(MLMG::BottomSolver::smoother);
   } else if (bottom_solver_type == "hypre") { 
      solver.setBottomSolver(MLMG::BottomSolver::hypre);
   }

    solver.setMaxIter(mg_max_iter);
    solver.setMaxFmgIter(mg_max_fmg_iter);
    solver.setVerbose(mg_verbose);
    solver.setCGVerbose(mg_cg_verbose);
    solver.setCGMaxIter(mg_cg_maxiter);

    // This ensures phi ghost cells are correctly filled when returned from solver
    solver.setFinalFillBC(true);

    //
    // Finally, solve the system
    //
    solver.solve(GetVecOfPtrs(this_phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol);
    solver.getFluxes(GetVecOfPtrs(fluxes), MLMG::Location::CellCenter);

    this_phi[lev]->FillBoundary(geom[lev].periodicity());
}

//
// Computes bcoeff = 1/ro at the faces of the scalar cells
//
void incflo::incflo_compute_bcoeff_ppe()
{
	BL_PROFILE("incflo::incflo_compute_bcoeff_ppe");

	// Directions
	int xdir = 1;
	int ydir = 2;
	int zdir = 3;

    for(int lev = 0; lev < nlev; lev++)
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
        bcoeff[lev][0]->FillBoundary(geom[lev].periodicity());
        bcoeff[lev][1]->FillBoundary(geom[lev].periodicity());
        bcoeff[lev][2]->FillBoundary(geom[lev].periodicity());
    }
}
