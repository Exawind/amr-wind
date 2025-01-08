#include <AMReX_BC_TYPES.H>
#include <memory>
#include "amr-wind/incflo.H"
#include "amr-wind/core/MLMGOptions.H"
#include "amr-wind/utilities/console_io.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/projection/nodal_projection_ops.H"
#include "hydro_utils.H"

using namespace amrex;

void amr_wind::nodal_projection::set_inflow_velocity(
    amr_wind::PhysicsMgr& phy_mgr,
    amr_wind::Field& vel_fld,
    int lev,
    amrex::Real time,
    MultiFab& vel_mfab,
    int nghost)
{
    vel_fld.set_inflow(lev, time, vel_mfab, nghost);

    // TODO fix hack for ABL
    if (phy_mgr.contains("ABL")) {
        auto& abl = phy_mgr.get<amr_wind::ABL>();
        const auto& bndry_plane = abl.bndry_plane();
        bndry_plane.populate_data(lev, time, vel_fld, vel_mfab);
        abl.abl_mpl().set_velocity(lev, time, vel_fld, vel_mfab);
    }
    if (phy_mgr.contains("OceanWaves")) {
        auto& ow = phy_mgr.get<amr_wind::ocean_waves::OceanWaves>();
        ow.ow_bndry().set_velocity(lev, time, vel_fld, vel_mfab);
    }
}

Array<amrex::LinOpBCType, AMREX_SPACEDIM>
amr_wind::nodal_projection::get_projection_bc(
    Orientation::Side side,
    amr_wind::Field& pressure,
    const Array<int, AMREX_SPACEDIM>& is_periodic)
{
    const auto& bctype = pressure.bc_type();
    Array<LinOpBCType, AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (is_periodic[dir] == 1) {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc = bctype[Orientation(dir, side)];
            switch (bc) {
            case BC::pressure_outflow: {
                r[dir] = LinOpBCType::Dirichlet;
                break;
            }
            case BC::mass_inflow_outflow:
            case BC::mass_inflow:
            case BC::wave_generation: {
                r[dir] = LinOpBCType::inflow;
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

void amr_wind::nodal_projection::apply_dirichlet_vel(
    amrex::MultiFab& mf_velocity, amrex::iMultiFab& mf_iblank)
{
    const auto& vel = mf_velocity.arrays();
    const auto& iblank = mf_iblank.const_arrays();

    amrex::ParallelFor(
        mf_velocity, mf_velocity.n_grow, mf_velocity.n_comp,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
            // Pure solid-body points
            if (iblank[nbx](i, j, k) == 0) {
                // Set velocity to 0 for now
                vel[nbx](i, j, k, n) = 0.0;
            }
        });
}

void amr_wind::nodal_projection::enforce_inout_solvability(
    amr_wind::Field& velocity,
    const Vector<Geometry>& geom,
    const int num_levels)
{
    BCRec const* bc_type = velocity.bcrec().data();
    Vector<Array<MultiFab*, AMREX_SPACEDIM>> vel_vec(num_levels);

    for (int lev = 0; lev < num_levels; ++lev) {
        vel_vec[lev][0] = new MultiFab(velocity(lev), amrex::make_alias, 0, 1);
        vel_vec[lev][1] = new MultiFab(velocity(lev), amrex::make_alias, 1, 1);
        vel_vec[lev][2] = new MultiFab(velocity(lev), amrex::make_alias, 2, 1);
    }

    HydroUtils::enforceInOutSolvability(vel_vec, bc_type, geom, true);
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
        (time > 0.0 and m_time.delta_t() < 0.1 * m_time.delta_t_nm1());

    if (m_verbose > 2) {
        if (proj_for_small_dt) {
            PrintMaxValues("before projection (small dt mod)");
        } else {
            PrintMaxValues("before projection");
        }
    }

    const bool is_anelastic = m_sim.is_anelastic();
    const bool variable_density =
        (!m_sim.pde_manager().constant_density() ||
         m_sim.physics_manager().contains("MultiPhase"));

    bool mesh_mapping = m_sim.has_mesh_mapping();

    auto& grad_p = m_repo.get_field("gp");
    auto& pressure = m_repo.get_field("p");
    auto& velocity = icns().fields().field;
    auto& velocity_old = icns().fields().field.state(amr_wind::FieldState::Old);
    amr_wind::Field const* mesh_fac =
        mesh_mapping
            ? &(m_repo.get_mesh_mapping_field(amr_wind::FieldLoc::CELL))
            : nullptr;
    amr_wind::Field const* mesh_detJ =
        mesh_mapping
            ? &(m_repo.get_mesh_mapping_det_j(amr_wind::FieldLoc::CELL))
            : nullptr;
    const auto* ref_density =
        is_anelastic ? &(m_repo.get_field("reference_density")) : nullptr;

    // TODO: Mesh mapping doesn't work with immersed boundaries
    // Do the pre pressure correction work -- this applies to IB only
    for (auto& pp : m_sim.physics()) {
        pp->pre_pressure_correction_work();
    }

    // ensure velocity is in stretched mesh space
    if (velocity.in_uniform_space() && mesh_mapping) {
        velocity.to_stretched_space();
    }

    // Add the ( grad p /ro ) back to u* (note the +dt)
    // Also account for mesh mapping in ( grad p /ro ) ->  1/fac * grad(p) *
    // dt/rho
    if (!incremental) {
        for (int lev = 0; lev <= finest_level; lev++) {

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(velocity(lev), TilingIfNotGPU()); mfi.isValid();
                 ++mfi) {
                Box const& bx = mfi.tilebox();
                Array4<Real> const& u = velocity(lev).array(mfi);
                Array4<Real const> const& rho = density[lev]->const_array(mfi);
                Array4<Real const> const& gp = grad_p(lev).const_array(mfi);
                amrex::Array4<amrex::Real const> fac =
                    mesh_mapping ? ((*mesh_fac)(lev).const_array(mfi))
                                 : amrex::Array4<amrex::Real const>();

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        Real soverrho = scaling_factor / rho(i, j, k);
                        amrex::Real fac_x =
                            mesh_mapping ? (fac(i, j, k, 0)) : 1.0;
                        amrex::Real fac_y =
                            mesh_mapping ? (fac(i, j, k, 1)) : 1.0;
                        amrex::Real fac_z =
                            mesh_mapping ? (fac(i, j, k, 2)) : 1.0;

                        u(i, j, k, 0) += 1 / fac_x * gp(i, j, k, 0) * soverrho;
                        u(i, j, k, 1) += 1 / fac_y * gp(i, j, k, 1) * soverrho;
                        u(i, j, k, 2) += 1 / fac_z * gp(i, j, k, 2) * soverrho;
                    });
            }
        }
    }

    // ensure velocity is in stretched mesh space
    if (velocity_old.in_uniform_space() && mesh_mapping) {
        velocity_old.to_stretched_space();
    }

    // Define "vel" to be U^* - U^n rather than U^*
    if (proj_for_small_dt || incremental) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Subtract(
                velocity(lev), velocity_old(lev), 0, 0, AMREX_SPACEDIM, 0);
        }
    }

    // scale U^* to accommodate for mesh mapping -> U^bar = J/fac * U
    if (mesh_mapping) {
        velocity.to_uniform_space();
    }

    // Create sigma while accounting for mesh mapping
    // sigma = 1/(fac^2)*J * dt/rho
    Vector<amrex::MultiFab> sigma(finest_level + 1);
    if (variable_density || mesh_mapping) {
        int ncomp = mesh_mapping ? AMREX_SPACEDIM : 1;
        for (int lev = 0; lev <= finest_level; ++lev) {
            sigma[lev].define(
                grids[lev], dmap[lev], ncomp, 0, MFInfo(), Factory(lev));
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(sigma[lev], TilingIfNotGPU()); mfi.isValid();
                 ++mfi) {
                Box const& bx = mfi.tilebox();
                Array4<Real> const& sig = sigma[lev].array(mfi);
                Array4<Real const> const& rho = density[lev]->const_array(mfi);
                amrex::Array4<amrex::Real const> fac =
                    mesh_mapping ? ((*mesh_fac)(lev).const_array(mfi))
                                 : amrex::Array4<amrex::Real const>();
                amrex::Array4<amrex::Real const> detJ =
                    mesh_mapping ? ((*mesh_detJ)(lev).const_array(mfi))
                                 : amrex::Array4<amrex::Real const>();
                const auto& ref_rho = is_anelastic
                                          ? (*ref_density)(lev).const_array(mfi)
                                          : amrex::Array4<amrex::Real>();

                amrex::ParallelFor(
                    bx, ncomp,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                        amrex::Real fac_cc =
                            mesh_mapping ? (fac(i, j, k, n)) : 1.0;
                        amrex::Real det_j =
                            mesh_mapping ? (detJ(i, j, k)) : 1.0;
                        sig(i, j, k, n) = std::pow(fac_cc, -2.) * det_j *
                                          scaling_factor / rho(i, j, k);
                        if (is_anelastic) {
                            sig(i, j, k, n) *= ref_rho(i, j, k);
                        }
                    });
            }
        }
    }

    // Perform projection
    std::unique_ptr<Hydro::NodalProjector> nodal_projector;

    auto bclo = amr_wind::nodal_projection::get_projection_bc(
        Orientation::low, pressure, m_sim.mesh().Geom()[0].isPeriodic());
    auto bchi = amr_wind::nodal_projection::get_projection_bc(
        Orientation::high, pressure, m_sim.mesh().Geom()[0].isPeriodic());

    Vector<MultiFab*> vel;
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel.push_back(&(velocity(lev)));
        vel[lev]->setBndry(0.0);
        if (!proj_for_small_dt and !incremental) {
            amr_wind::nodal_projection::set_inflow_velocity(
                m_sim.physics_manager(), velocity, lev, time, *vel[lev], 1);

            // fill periodic boundaries to avoid corner cell issues
            vel[lev]->FillBoundary(geom[lev].periodicity());
        }
    }

    // Need to apply custom Neumann funcs for inflow-outflow BC
    // after setting the inflow vels above
    // and then enforce solvability by matching outflow to inflow.
    if (!proj_for_small_dt and !incremental and velocity.has_inout_bndry()) {
        velocity.apply_bc_funcs(amr_wind::FieldState::New);

        amr_wind::nodal_projection::enforce_inout_solvability(
            velocity, m_repo.mesh().Geom(), m_repo.num_active_levels());
    }

    if (is_anelastic) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (int idim = 0; idim < velocity.num_comp(); ++idim) {
                amrex::Multiply(
                    velocity(lev), (*ref_density)(lev), 0, idim,
                    density[lev]->nComp(), 0);
            }
        }
    }

    amr_wind::MLMGOptions options("nodal_proj");

    if (variable_density || mesh_mapping) {
        nodal_projector = std::make_unique<Hydro::NodalProjector>(
            vel, GetVecOfConstPtrs(sigma), Geom(0, finest_level),
            options.lpinfo());
    } else {
        amrex::Real rho_0 = 1.0;
        amrex::ParmParse pp("incflo");
        pp.query("density", rho_0);

        nodal_projector = std::make_unique<Hydro::NodalProjector>(
            vel, scaling_factor / rho_0, Geom(0, finest_level),
            options.lpinfo());
    }

    // Set MLMG and NodalProjector options
    options(*nodal_projector);
    nodal_projector->setDomainBC(bclo, bchi);

    bool has_ib = m_sim.physics_manager().contains("IB");
    if (has_ib) {
        auto div_vel_rhs =
            sim().repo().create_scratch_field(1, 0, amr_wind::FieldLoc::NODE);
        nodal_projector->computeRHS(div_vel_rhs->vec_ptrs(), vel, {}, {});
        // Mask the righ-hand side of the Poisson solve for the nodes inside the
        // body
        const auto& imask_node = repo().get_int_field("mask_node");
        for (int lev = 0; lev <= finest_level; ++lev) {
            amrex::MultiFab::Multiply(
                *div_vel_rhs->vec_ptrs()[lev],
                amrex::ToMultiFab(imask_node(lev)), 0, 0, 1, 0);
        }
        nodal_projector->setCustomRHS(div_vel_rhs->vec_const_ptrs());
    }

    if (m_sim.has_overset()) {
        // Setup masking for overset simulations
        auto& linop = nodal_projector->getLinOp();
        const auto& imask_node = repo().get_int_field("mask_node");
        for (int lev = 0; lev <= finest_level; ++lev) {
            linop.setOversetMask(lev, imask_node(lev));
        }

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

    if (is_anelastic) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (int idim = 0; idim < velocity.num_comp(); ++idim) {
                amrex::Divide(
                    velocity(lev), (*ref_density)(lev), 0, idim,
                    density[lev]->nComp(), 0);
            }
        }
    }

    // scale U^* back to -> U = fac/J * U^bar
    if (mesh_mapping) {
        velocity.to_stretched_space();
    }

    // Define "vel" to be U^{n+1} rather than (U^{n+1}-U^n)
    if (proj_for_small_dt || incremental) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Add(
                velocity(lev), velocity_old(lev), 0, 0, AMREX_SPACEDIM, 0);
        }
    }

    // Get phi and fluxes
    auto phi = nodal_projector->getPhi();
    auto gradphi = nodal_projector->getGradPhi();

    for (int lev = 0; lev <= finest_level; lev++) {

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
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

    // Determine if reference pressure should be added back
    if (m_reconstruct_true_pressure && time != 0.0) {
        const auto& p0 = m_repo.get_field("reference_pressure");
        for (int lev = 0; lev <= finest_level; lev++) {
            amrex::MultiFab::Add(
                pressure(lev), p0(lev), 0, 0, 1, p0.num_grow()[0]);
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
