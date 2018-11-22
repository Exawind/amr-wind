#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

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
// except in the initial projection when
//
//     new p  = old p + phi     (nstep has its initial value -1)
//
// Note: scaling_factor equals dt except when called during initial projection, when it is 1.0
//
void incflo::ApplyProjection(Real time, Real scaling_factor)
{
	BL_PROFILE("incflo::ApplyProjection");

    // Print info about predictor step
    if(incflo_verbose > 1)
    {
        amrex::Print() << "Before projection:" << std::endl;
        PrintMaxValues(time);
    }

    // Add the (1/rho gradp) back to ustar (note the +dt)
    if(nstep >= 0)
    {
        for(int lev = 0; lev <= finest_level; lev++)
        {
            // Convert velocities to momenta
            for(int dir = 0; dir < 3; dir++)
            {
                MultiFab::Multiply(*vel[lev], *ro[lev], 0, dir, 1, vel[lev]->nGrow());
            }

            MultiFab::Saxpy(*vel[lev], scaling_factor, *gp[lev], 0, 0, 3, vel[lev]->nGrow());

            // Convert momenta back to velocities
            for(int dir = 0; dir < 3; dir++)
            {
                MultiFab::Divide(*vel[lev], *ro[lev], 0, dir, 1, vel[lev]->nGrow());
            }
        }
    }

    // Make sure div(u) is updated before scaling by dt (to be used as RHS for Poisson equation)
    ComputeDivU(time);
    for(int lev = 0; lev <= finest_level; lev++)
    {
        divu[lev]->mult(1.0 / scaling_factor, divu[lev]->nGrow());
    }

    // Declare, resize, reset and initialize MultiFabs to hold the solution of the Poisson solve
	Vector<std::unique_ptr<MultiFab>> phi;
	Vector<std::unique_ptr<MultiFab>> fluxes;
    phi.resize(finest_level + 1);
    fluxes.resize(finest_level + 1);
    for(int lev = 0; lev <= finest_level; lev++)
    {
        const BoxArray & nd_grids = amrex::convert(grids[lev], IntVect{1,1,1});
        phi[lev].reset(new MultiFab(nd_grids, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));
        phi[lev]->setVal(0.0);
        fluxes[lev].reset(new MultiFab(vel[lev]->boxArray(),
                                       vel[lev]->DistributionMap(),
                                       vel[lev]->nComp(), 1,
                                       MFInfo(), *ebfactory[lev]));
        fluxes[lev]->setVal(1.0e200);
    }

    //
    // Solve Poisson Equation:
    //
    //                  div( 1/rho * grad(phi) ) = divu / dt
    //
    // Also outputs minus grad(phi) / rho into "fluxes"
    //
	poisson_equation->solve(phi, fluxes, ro, divu);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // Multiply fluxes by dt
        fluxes[lev]->mult(scaling_factor, fluxes[lev]->nGrow());

        // Now we correct the velocity with MINUS dt * (1/rho) * grad(phi),
        MultiFab::Add(*vel[lev], *fluxes[lev], 0, 0, 3, 0);

        // Multiply by rho and divide by (-dt) to get grad(phi)
        fluxes[lev]->mult(-1.0 / scaling_factor, fluxes[lev]->nGrow());
        for(int dir = 0; dir < 3; dir++)
        {
            MultiFab::Multiply(*fluxes[lev], (*ro[lev]), 0, dir, 1, fluxes[lev]->nGrow());
        }

        if(nstep >= 0)
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

    // Print info about predictor step
    if(incflo_verbose > 1)
    {
        amrex::Print() << "After projection: " << std::endl; 
        PrintMaxValues(time);
    }
}
