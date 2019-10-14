#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

#include <incflo.H>
#include <NodalProjection.H>

//
// Computes the following decomposition:
//
//    u + dt grad(phi) / ro = u*,     where div(u) = 0
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

    if(incflo_verbose > 2)
    {
        amrex::Print() << "Before projection:" << std::endl;
        PrintMaxValues(time);
    }

    // Add the ( grad p /ro ) back to u* (note the +dt)
    if(nstep >= 0)
    {
        for(int lev = 0; lev <= finest_level; lev++)
        {
            // Convert velocities to momenta
            for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
            {
                MultiFab::Multiply(*vel[lev], *density[lev], 0, dir, 1, vel[lev]->nGrow());
            }

            MultiFab::Saxpy(*vel[lev], scaling_factor, *gp[lev], 0, 0, AMREX_SPACEDIM, vel[lev]->nGrow());

            // Convert momenta back to velocities
            for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
            {
                MultiFab::Divide(*vel[lev], *density[lev], 0, dir, 1, vel[lev]->nGrow());
            }
        }
    }

    nodal_projector -> project(vel, density, time, scaling_factor);

    // Get phi and fluxes
    Vector< const amrex::MultiFab* >  phi(nlev);
    Vector< const amrex::MultiFab* >  gradphi(nlev);

    phi     = nodal_projector -> getPhi();
    gradphi = nodal_projector -> getGradPhi();

    // // Make sure div(u) is up to date
    // ComputeDivU(time);

    // // Declare, resize, reset and initialize MultiFabs to hold the solution of the Poisson solve
    // Vector<std::unique_ptr<MultiFab>> phi;
    // Vector<std::unique_ptr<MultiFab>> fluxes;

    // phi.resize(finest_level + 1);
    // fluxes.resize(finest_level + 1);

    // for(int lev = 0; lev <= finest_level; lev++)
    // {
    //     const BoxArray & nd_grids = amrex::convert(grids[lev], IntVect{1,1,1});
    //     phi[lev].reset(new MultiFab(nd_grids, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));
    //     phi[lev]->setVal(0.0);
    //     fluxes[lev].reset(new MultiFab(vel[lev]->boxArray(),
    //                                    vel[lev]->DistributionMap(),
    //                                    vel[lev]->nComp(), 1,
    //                                    MFInfo(), *ebfactory[lev]));
    //     fluxes[lev]->setVal(1.0e200);
    // }

    // //
    // // Solve Poisson Equation:
    // //
    // //                  div( 1/rho * grad(phi) ) = divu
    // //
    // // Note that
    // //      p = phi / dt
    // // for all steps except the initial projection, when we add phi / dt to p instead.
    // //
    // // Also outputs minus grad(phi) / rho into "fluxes"
    // //
    // poisson_equation->solve(phi, fluxes, density, divu);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        if(nstep >= 0)
        {
            // p := phi
            MultiFab::Copy(*p[lev], *phi[lev], 0, 0, 1, phi[lev]->nGrow());
            MultiFab::Copy(*gp[lev], *gradphi[lev], 0, 0, AMREX_SPACEDIM, gradphi[lev]->nGrow());
        }
        else
        {
            // p := p + phi
            MultiFab::Add(*p[lev], *phi[lev], 0, 0, 1, phi[lev]->nGrow());
            MultiFab::Add(*gp[lev], *gradphi[lev], 0, 0, AMREX_SPACEDIM, gradphi[lev]->nGrow());
        }
    }

    AverageDown();

    if(incflo_verbose > 2)
    {
        amrex::Print() << "After projection: " << std::endl;
        PrintMaxValues(time);
    }
}
