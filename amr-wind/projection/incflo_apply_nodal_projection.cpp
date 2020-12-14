#include <AMReX_BC_TYPES.H>
#include "amr-wind/incflo.H"
#include "amr-wind/core/MLMGOptions.H"
#include "amr-wind/utilities/console_io.H"
#include "amr-wind/core/field_ops.H"

using namespace amrex;

Array<amrex::LinOpBCType, AMREX_SPACEDIM>
incflo::get_projection_bc(Orientation::Side side) const noexcept
{
    auto& bctype = pressure().bc_type();
    Array<LinOpBCType, AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (geom[0].isPeriodic(dir)) {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc = bctype[Orientation(dir, side)];
            switch (bc) {
            case BC::pressure_inflow:
            case BC::pressure_outflow: {
                r[dir] = LinOpBCType::Dirichlet;
                break;
            }
            default:
                r[dir] = LinOpBCType::Neumann;
                break;
            };
        }
    }
    return r;
}

/** Perform nodal projection
 *
 *  Computes the following decomposition:
 *
 *  \f{align}
 *    u + \frac{\Delta t}{\rho} \nabla \phi &= u^{*} & \nabla \cdot u = 0
 *  \f}
 *
 *  where \f$u^{*}\f$ is a non-div-free velocity field, stored
 *  by components in u, v, and w. The resulting div-free
 *  velocity field, u, overwrites the value of u* in u, v, and w.
 *
 *  \f$\phi\f$ is an auxiliary function related to the pressure \f$p\f$ by the
 * relation:
 *
 *  \f[ p^{\mathrm{new}} = \phi \f]
 *
 *  except in the initial projection when
 *
 *  \f[ p^{\mathrm{new}} = p^{\mathrm{old}} + \phi \f]
 *
 *  Velocity is updated using the following relation
 *
 *  \f{align}
 *   u^{**} &= u^{*} + \frac{\Delta t}{\rho} \nabla p^{\mathrm{old}} \\
 *   \nabla \cdot \left( \frac{\Delta t}{\rho} \nabla \phi \right) &= \nabla
 * \cdot u^{**} \\ u^{\mathrm{new}} &= u^{**} - \frac{\Delta t}{\rho} \nabla
 * p^\mathrm{new} \f}
 *
 *  Notes:
 *  - `scaling_factor` equals \f$\Delta t\f$ except when called during initial
 *     projection, when it is 1.0
 *
 *  - \f$\rho\f$ in the above expressions is either at state `n+1` or `n+1/2`
 *    depending on whether this method was called from incflo::ApplyPredictor or
 *    incflo::ApplyCorrector.
 *
 *  - If `incremental == true`, then the pressure term is not added to
 *    \f$u^{**}\f$ and the update is in delta-form.
 *
 *  Please consult [AMReX Linear
 *  Solvers](https://amrex-codes.github.io/amrex/docs_html/LinearSolvers.html#nodal-projection)
 *  documentation for more information on the nodal projection operator.
 *
 *  \param density Vector of multifabs containing the density field at all
 * levels \param time Current time \param scaling_factor Scaling factor
 * \f$\Delta t\f$ \param incremental Flag indicating if it is an incremental
 * projection.
 *
 */
void incflo::ApplyProjection(
    Vector<MultiFab const*> density,
    Real time,
    Real scaling_factor,
    bool incremental)
{
    BL_PROFILE("amr-wind::incflo::ApplyProjection");

    // If we have dropped the dt substantially for whatever reason,
    // use a different form of the approximate projection that
    // projects (U^*-U^n + dt Gp) rather than (U^* + dt Gp)

    bool proj_for_small_dt =
        (time > 0.0 and m_time.deltaT() < 0.1 * m_time.deltaTNm1());

    if (m_verbose > 2) {
        if (proj_for_small_dt) {
            PrintMaxValues("before projection (small dt mod)");
        } else {
            PrintMaxValues("before projection");
        }
    }

    bool variable_density =
        (!m_sim.pde_manager().constant_density() ||
         m_sim.physics_manager().contains("MultiPhase"));

    auto& grad_p = m_repo.get_field("gp");
    auto& pressure = m_repo.get_field("p");
    auto& velocity = icns().fields().field;

    // Add the ( grad p /ro ) back to u* (note the +dt)
    if (!incremental) {
        for (int lev = 0; lev <= finest_level; lev++) {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(velocity(lev), TilingIfNotGPU()); mfi.isValid();
                 ++mfi) {
                Box const& bx = mfi.tilebox();
                Array4<Real> const& u = velocity(lev).array(mfi);
                Array4<Real const> const& rho = density[lev]->const_array(mfi);
                Array4<Real const> const& gp = grad_p(lev).const_array(mfi);
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        Real soverrho = scaling_factor / rho(i, j, k);
                        u(i, j, k, 0) += gp(i, j, k, 0) * soverrho;
                        u(i, j, k, 1) += gp(i, j, k, 1) * soverrho;
                        u(i, j, k, 2) += gp(i, j, k, 2) * soverrho;
                    });
            }
        }
    }

    // Define "vel" to be U^* - U^n rather than U^*
    if (proj_for_small_dt || incremental) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Subtract(
                velocity(lev), velocity.state(amr_wind::FieldState::Old)(lev),
                0, 0, AMREX_SPACEDIM, 0);
        }
    }

    // Create sigma
    Vector<amrex::MultiFab> sigma(finest_level + 1);
    if (variable_density) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            sigma[lev].define(
                grids[lev], dmap[lev], 1, 0, MFInfo(), Factory(lev));
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(sigma[lev], TilingIfNotGPU()); mfi.isValid();
                 ++mfi) {
                Box const& bx = mfi.tilebox();
                Array4<Real> const& sig = sigma[lev].array(mfi);
                Array4<Real const> const& rho = density[lev]->const_array(mfi);
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        sig(i, j, k) = scaling_factor / rho(i, j, k);
                    });
            }
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
            //  velocity.fillphysbc(lev, time, *vel[lev], 1);
        }
    }

    amr_wind::MLMGOptions options("nodal_proj");

    if (variable_density) {
        nodal_projector.reset(new NodalProjector(
            vel, GetVecOfConstPtrs(sigma), Geom(0, finest_level),
            options.lpinfo()));
    } else {

        amrex::Real rho_0 = 1.0;
        amrex::ParmParse pp("incflo");
        pp.query("density", rho_0);

        nodal_projector.reset(new NodalProjector(
            vel, scaling_factor / rho_0, Geom(0, finest_level),
            options.lpinfo()));
    }

    // Set MLMG and NodalProjector options
    options(*nodal_projector);
    nodal_projector->setDomainBC(bclo, bchi);

    // Setup masking for overset simulations
    if (sim().has_overset()) {
        auto& linop = nodal_projector->getLinOp();
        auto& imask_node = repo().get_int_field("mask_node");
        for (int lev = 0; lev <= finest_level; ++lev) {
            linop.setOversetMask(lev, imask_node(lev));
        }
    }

    if (m_sim.has_overset()) {
        auto phif = m_repo.create_scratch_field(1, 1, amr_wind::FieldLoc::NODE);
        if (incremental) {
            for (int lev = 0; lev <= finestLevel(); ++lev) {
                (*phif)(lev).setVal(0.0);
            }
        } else {
            amr_wind::field_ops::copy(*phif, pressure, 0, 0, 1, 1);
        }

        nodal_projector->project(
            phif->vec_ptrs(), options.rel_tol, options.abs_tol);
    } else {
        nodal_projector->project(options.rel_tol, options.abs_tol);
    }
    amr_wind::io::print_mlmg_info(
        "Nodal_projection", nodal_projector->getMLMG());

    // Define "vel" to be U^{n+1} rather than (U^{n+1}-U^n)
    if (proj_for_small_dt || incremental) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Add(
                velocity(lev), velocity.state(amr_wind::FieldState::Old)(lev),
                0, 0, AMREX_SPACEDIM, 0);
        }
    }

    // Get phi and fluxes
    auto phi = nodal_projector->getPhi();
    auto gradphi = nodal_projector->getGradPhi();

    for (int lev = 0; lev <= finest_level; lev++) {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(grad_p(lev), TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& tbx = mfi.tilebox();
            Box const& nbx = mfi.nodaltilebox();
            Array4<Real> const& gp_lev = grad_p(lev).array(mfi);
            Array4<Real> const& p_lev = pressure(lev).array(mfi);
            Array4<Real const> const& gp_proj = gradphi[lev]->const_array(mfi);
            Array4<Real const> const& p_proj = phi[lev]->const_array(mfi);
            if (incremental) {
                amrex::ParallelFor(
                    tbx, AMREX_SPACEDIM,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                        gp_lev(i, j, k, n) += gp_proj(i, j, k, n);
                    });
                amrex::ParallelFor(
                    nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        p_lev(i, j, k) += p_proj(i, j, k);
                    });
            } else {
                amrex::ParallelFor(
                    tbx, AMREX_SPACEDIM,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                        gp_lev(i, j, k, n) = gp_proj(i, j, k, n);
                    });
                amrex::ParallelFor(
                    nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        p_lev(i, j, k) = p_proj(i, j, k);
                    });
            }
        }
    }

    for (int lev = finest_level - 1; lev >= 0; --lev) {
        amrex::average_down(
            grad_p(lev + 1), grad_p(lev), 0, AMREX_SPACEDIM, refRatio(lev));
    }

    velocity.fillpatch(m_time.new_time());
    if (m_verbose > 2) {
        if (proj_for_small_dt) {
            PrintMaxValues("after projection (small dt mod)");
        } else {
            PrintMaxValues("after projection");
        }
    }
}
