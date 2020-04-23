#include <AMReX_BC_TYPES.H>
#include <incflo.H>
#include <MLMGOptions.H>
#include "console_io.H"

using namespace amrex;

Array<amrex::LinOpBCType,AMREX_SPACEDIM>
incflo::get_projection_bc (Orientation::Side side) const noexcept
{
    auto& bctype = pressure().bc_type();
    Array<LinOpBCType,AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (geom[0].isPeriodic(dir)) {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc = bctype[Orientation(dir,side)];
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
            case BC::wall_model:
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
void incflo::ApplyProjection (Vector<MultiFab const*> density,
                              Real time, Real scaling_factor, bool incremental)
{
    BL_PROFILE("amr-wind::incflo::ApplyProjection")

    // If we have dropped the dt substantially for whatever reason,
    // use a different form of the approximate projection that
    // projects (U^*-U^n + dt Gp) rather than (U^* + dt Gp)

    bool proj_for_small_dt = (time > 0.0 and m_time.deltaT() < 0.1 * m_time.deltaTNm1());

    if (m_verbose > 2)
    {
        if (proj_for_small_dt) {
            amrex::Print() << "Before projection (with small dt modification):" << std::endl;
        } else {
            amrex::Print() << "Before projection:" << std::endl;
        }
        PrintMaxValues(time);
    }

    auto& grad_p = m_repo.get_field("gp");
    auto& pressure = m_repo.get_field("p");
    auto& velocity = icns().fields().field;

    // Add the ( grad p /ro ) back to u* (note the +dt)
    if (!incremental)
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(velocity(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real> const& u = velocity(lev).array(mfi);
                Array4<Real const> const& rho = density[lev]->const_array(mfi);
                Array4<Real const> const& gp = grad_p(lev).const_array(mfi);
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
            MultiFab::Subtract(velocity(lev),
                               velocity.state(amr_wind::FieldState::Old)(lev),
                               0, 0, AMREX_SPACEDIM, 0);
        }
    }

    // Create sigma
    Vector<amrex::MultiFab> sigma(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev )
    {
        sigma[lev].define(grids[lev], dmap[lev], 1, 0, MFInfo(), Factory(lev));
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(sigma[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real> const& sig = sigma[lev].array(mfi);
            Array4<Real const> const& rho = density[lev]->const_array(mfi);
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
        vel.push_back(&(velocity(lev)));
        vel[lev]->setBndry(0.0);
        if (!proj_for_small_dt and !incremental) {
            set_inflow_velocity(lev, time, *vel[lev], 1);
        }
    }

    amr_wind::MLMGOptions options("nodal_proj");
    nodal_projector.reset(new NodalProjector(vel, GetVecOfConstPtrs(sigma),
                                             Geom(0,finest_level), LPInfo()));
    nodal_projector->setDomainBC(bclo, bchi);
    nodal_projector->setVerbose(options.verbose);
    nodal_projector->project(options.rel_tol, options.abs_tol);
    amr_wind::io::print_mlmg_info("Nodal_projection", nodal_projector->getMLMG());

    // Define "vel" to be U^{n+1} rather than (U^{n+1}-U^n)
    if (proj_for_small_dt || incremental)
    {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Add(velocity(lev),
                          velocity.state(amr_wind::FieldState::Old)(lev), 0, 0, AMREX_SPACEDIM, 0);
        }
    }

    // Get phi and fluxes
    auto phi = nodal_projector->getPhi();
    auto gradphi = nodal_projector->getGradPhi();

    for(int lev = 0; lev <= finest_level; lev++)
    {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(grad_p(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& tbx = mfi.tilebox();
            Box const& nbx = mfi.nodaltilebox();
            Array4<Real> const& gp_lev = grad_p(lev).array(mfi);
            Array4<Real> const& p_lev = pressure(lev).array(mfi);
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
        amrex::average_down(grad_p(lev+1), grad_p(lev),
                            0, AMREX_SPACEDIM, refRatio(lev));
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
