#include <AMReX_BC_TYPES.H>
#include <incflo.H>

using namespace amrex;

Array<amrex::LinOpBCType,AMREX_SPACEDIM>
incflo::get_projection_bc (Orientation::Side side) const noexcept
{
    Array<LinOpBCType,AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (geom[0].isPeriodic(dir)) {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc = m_bc_type[Orientation(dir,side)];
            switch (bc)
            {
            case BC::pressure_inflow:
            case BC::pressure_outflow:
            {
                r[dir] = LinOpBCType::Dirichlet;
                break;
            }
            case BC::mass_inflow:
            case BC::slip_wall:
            case BC::no_slip_wall:
            {
                r[dir] = LinOpBCType::Neumann;
                break;
            }
            default:
                amrex::Abort("get_projection_bc: undefined BC type");
            };
        }
    }
    return r;
}

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
void incflo::ApplyProjection (Real time, Real scaling_factor, bool incremental)
{
    BL_PROFILE("incflo::ApplyProjection");

    // If we have dropped the dt substantially for whatever reason,
    // use a different form of the approximate projection that
    // projects (U^*-U^n + dt Gp) rather than (U^* + dt Gp)

    bool proj_for_small_dt = (time > 0.0 and m_dt < 0.1 * m_prev_dt);

    if (m_verbose > 2)
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
        for (int lev = 0; lev <= finest_level; lev++)
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
    if (proj_for_small_dt || incremental)
    {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Subtract(m_leveldata[lev]->velocity,
                               m_leveldata[lev]->velocity_o, 0, 0, AMREX_SPACEDIM, 0);
        }
    }

    // Create sigma
    Vector<amrex::MultiFab> sigma(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev )
    {
        auto const& ld = *m_leveldata[lev];
        sigma[lev].define(grids[lev], dmap[lev], 1, 0, MFInfo(), *m_factory[lev]);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(sigma[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real> const& sig = sigma[lev].array(mfi);
            Array4<Real const> const& rho = ld.density.const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                sig(i,j,k) = scaling_factor / rho(i,j,k);
            });
        }
    }

    // Perform projection
    std::unique_ptr<NodalProjector> nodal_projector;

    auto bclo = get_projection_bc(Orientation::low);
    auto bchi = get_projection_bc(Orientation::high);

    Vector<MultiFab*> vel;
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel.push_back(&(m_leveldata[lev]->velocity));
        vel[lev]->setBndry(0.0);
        if (!proj_for_small_dt and !incremental) {
            set_inflow_velocity(lev, time, *vel[lev], 1);
        }
    }

    nodal_projector.reset(new NodalProjector(vel, GetVecOfConstPtrs(sigma),
                                             Geom(0,finest_level), LPInfo()));
    nodal_projector->setDomainBC(bclo, bchi);
    nodal_projector->project();

    // Define "vel" to be U^{n+1} rather than (U^{n+1}-U^n)
    if (proj_for_small_dt || incremental)
    {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Add(m_leveldata[lev]->velocity,
                          m_leveldata[lev]->velocity_o, 0, 0, AMREX_SPACEDIM, 0);
        }
    }

    // Get phi and fluxes
    auto phi = nodal_projector->getPhi();
    auto gradphi = nodal_projector->getGradPhi();

    for(int lev = 0; lev <= finest_level; lev++)
    {
        auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(ld.gp,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& tbx = mfi.tilebox();
            Box const& nbx = mfi.nodaltilebox();
            Array4<Real> const& gp_lev = ld.gp.array(mfi);
            Array4<Real> const& p_lev = ld.p.array(mfi);
            Array4<Real const> const& gp_proj = gradphi[lev]->const_array(mfi);
            Array4<Real const> const& p_proj = phi[lev]->const_array(mfi);
            if (incremental) {
                amrex::ParallelFor(tbx, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    gp_lev(i,j,k,n) += gp_proj(i,j,k,n);
                });
                amrex::ParallelFor(nbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    p_lev (i,j,k) += p_proj(i,j,k);
                });
            } else {
                amrex::ParallelFor(tbx, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    gp_lev(i,j,k,n) = gp_proj(i,j,k,n);
                });
                amrex::ParallelFor(nbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    p_lev(i,j,k) = p_proj(i,j,k);
                });
            }
        }
    }

    for (int lev = finest_level-1; lev >= 0; --lev) {
#ifdef AMREX_USE_EB
        amrex::EB_average_down(m_leveldata[lev+1]->gp, m_leveldata[lev]->gp,
                               0, AMREX_SPACEDIM, refRatio(lev));
#else
        amrex::average_down(m_leveldata[lev+1]->gp, m_leveldata[lev]->gp,
                            0, AMREX_SPACEDIM, refRatio(lev));
#endif
    }

    if (m_verbose > 2)
    {
        if (proj_for_small_dt) {
            amrex::Print() << "After  projection (with small dt modification):" << std::endl;
        } else {
            amrex::Print() << "After  projection:" << std::endl;
        }
        PrintMaxValues(time);
    }
}
