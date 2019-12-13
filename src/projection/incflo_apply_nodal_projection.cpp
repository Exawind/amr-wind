#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

#include <incflo.H>

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
void incflo::ApplyProjection(Real time, Real scaling_factor, bool incremental)
{
    BL_PROFILE("incflo::ApplyProjection");

    // If we have dropped the dt substantially for whatever reason, use a different form of the approximate
    // projection that projects (U^*-U^n + dt Gp) rather than (U^* + dt Gp)

    bool proj_for_small_dt = (time > 0.0 && dt < 0.1 * prev_dt);

    if (incflo_verbose > 2)
    {
        if (proj_for_small_dt) {
            amrex::Print() << "Before projection (with small dt modification):" << std::endl;
        } else {
            amrex::Print() << "Before projection:" << std::endl;
        }
        PrintMaxValues(time);
    }

    // Add the ( grad p /ro ) back to u* (note the +dt)
    if (!incremental)
    {
        for(int lev = 0; lev <= finest_level; lev++)
        {
            auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real> const& u = ld.velocity.array(mfi);
                Array4<Real const> const& rho = ld.density.const_array(mfi);
                Array4<Real const> const& gp = ld.gp.const_array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real soverrho = scaling_factor / rho(i,j,k);
                    u(i,j,k,0) += gp(i,j,k,0) * soverrho;
                    u(i,j,k,1) += gp(i,j,k,1) * soverrho;
                    u(i,j,k,2) += gp(i,j,k,2) * soverrho;
                });
            }
        }
    }

    // Define "vel" to be U^* - U^n rather than U^*
#if 0
    // xxxxx todo
    if (proj_for_small_dt)
    {
       incflo_set_velocity_bcs(time, vel_o, extrap_dir_bcs);

       for(int lev = 0; lev <= finest_level; lev++)
          MultiFab::Saxpy(*vel[lev], -1.0, *vel_o[lev], 0, 0, AMREX_SPACEDIM, vel[lev]->nGrow());
    }
#endif

    // Create sigma
    Vector< std::unique_ptr< amrex::MultiFab > >  sigma(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev )
    {
        auto const& ld = *m_leveldata[lev];
        sigma[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *m_factory[lev]));
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*sigma[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real> const& sig = sigma[lev]->array(mfi);
            Array4<Real const> const& rho = ld.density.const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                sig(i,j,k) = scaling_factor / rho(i,j,k);
            });
        }
    }

    // Perform projection
    {
        if (!nodal_projector) {
            std::array<LinOpBCType,AMREX_SPACEDIM> bclo;
            std::array<LinOpBCType,AMREX_SPACEDIM> bchi;
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                if (geom[0].isPeriodic(dir)) {
                    bclo[dir] = LinOpBCType::Periodic;
                    bchi[dir] = LinOpBCType::Periodic;
                } else {
                    if (m_bc_type[Orientation(dir,Orientation::low)] == BC::pressure_inflow or
                        m_bc_type[Orientation(dir,Orientation::low)] == BC::pressure_outflow) {
                        bclo[dir] = LinOpBCType::Dirichlet;
                    } else {
                        bclo[dir] = LinOpBCType::Neumann;
                    }
                    if (m_bc_type[Orientation(dir,Orientation::high)] == BC::pressure_inflow or
                        m_bc_type[Orientation(dir,Orientation::high)] == BC::pressure_outflow) {
                        bchi[dir] = LinOpBCType::Dirichlet;
                    } else {
                        bchi[dir] = LinOpBCType::Neumann;
                    }
                }
            }
#ifdef AMREX_USE_EB
            if (!EBFactory(0).isAllRegular()) {
                Vector<EBFArrayBoxFactory const*> fact(finest_level+1);
                for (int lev = 0; lev <= finest_level; ++lev) {
                    fact[lev] = &(EBFactory(lev));
                }
                nodal_projector.reset(new NodalProjector(Geom(0,finest_level),
                                                         boxArray(0,finest_level),
                                                         DistributionMap(0,finest_level),
                                                         bclo, bchi, fact, LPInfo()));
            } else
#endif
            {
                nodal_projector.reset(new NodalProjector(Geom(0,finest_level),
                                                         boxArray(0,finest_level),
                                                         DistributionMap(0,finest_level),
                                                         bclo, bchi, LPInfo()));
            }
        }

        Vector<MultiFab*> vel;
        for (int lev = 0; lev <= finest_level; ++lev) {
            vel.push_back(&(m_leveldata[lev]->velocity));
            vel[lev]->setBndry(0.0);
            set_inflow_velocity(lev, time, *vel[lev], 1);
        }

        EB_set_covered(*vel[0], 0, 3, nghost, 0.0);
        EB_set_covered(*sigma[0], 0, 1, 0, 0.0);
        VisMF::Write(*vel[0], "b-vel");
        VisMF::Write(*sigma[0], "b-sig");
        amrex::Abort("xxxxx");

        nodal_projector->project(vel, GetVecOfConstPtrs(sigma));
    }

#if 0
    // xxxxx todo
    // Define "vel" to be U^{n+1} rather than (U^{n+1}-U^n)
    if (proj_for_small_dt)
    {
       for(int lev = 0; lev <= finest_level; lev++)
          MultiFab::Saxpy(*vel[lev], 1.0, *vel_o[lev], 0, 0, AMREX_SPACEDIM, vel[lev]->nGrow());
    }
#endif

    // Get phi and fluxes
    auto phi = nodal_projector->getPhi();
    auto gradphi = nodal_projector->getGradPhi();

    EB_set_covered(m_leveldata[0]->velocity, 0, 3, nghost, 0.0);
    VisMF::Write(m_leveldata[0]->velocity, "vel");
    VisMF::Write(*phi[0], "phi");
    VisMF::Write(*gradphi[0], "gradphi");

    amrex::Abort("xxxxx after project: so far so good");

    for(int lev = 0; lev <= finest_level; lev++)
    {
        if (incremental)
        {
            // p := p + phi
            MultiFab::Add(*p[lev], *phi[lev], 0, 0, 1, phi[lev]->nGrow());
            MultiFab::Add(*gp[lev], *gradphi[lev], 0, 0, AMREX_SPACEDIM, gradphi[lev]->nGrow());
        }
        else
        {
            // p := phi
            MultiFab::Copy(*p[lev], *phi[lev], 0, 0, 1, phi[lev]->nGrow());
            MultiFab::Copy(*gp[lev], *gradphi[lev], 0, 0, AMREX_SPACEDIM, gradphi[lev]->nGrow());
        }

    }

    AverageDown();

    if(incflo_verbose > 2)
    {
        if (proj_for_small_dt)
           amrex::Print() << "After  projection (with small dt modification):" << std::endl;
        else
           amrex::Print() << "After  projection:" << std::endl;
        PrintMaxValues(time);
    }
}
