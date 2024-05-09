#include "amr-wind/overset/OversetOps.H"
#include "amr-wind/core/field_ops.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFabUtil.H"
#include "amr-wind/core/MLMGOptions.H"
#include "amr-wind/projection/nodal_projection_ops.H"
#include <hydro_NodalProjector.H>

namespace amr_wind {

void OversetOps::initialize(CFDSim& sim) { m_sim_ptr = &sim; }

void OversetOps::pre_advance_actions()
{
    // Update pressure gradient using updated overset pressure field
    update_gradp();
}

/* ----------------------------------------------- */
/* PUBLIC FUNCTIONS ABOVE, PRIVATE FUNCTIONS BELOW */
/* ----------------------------------------------- */

void OversetOps::update_gradp()
{
    BL_PROFILE("amr-wind::OversetOps::update_gradp");

    auto density = (*m_sim_ptr).repo().get_field("density").vec_const_ptrs();
    auto scaling_factor = (*m_sim_ptr).time().deltaT();

    // Pressure and sigma are necessary to calculate the pressure gradient

    const bool is_anelastic = (*m_sim_ptr).is_anelastic();
    const bool variable_density =
        (!(*m_sim_ptr).pde_manager().constant_density() ||
         (*m_sim_ptr).physics_manager().contains("MultiPhase"));

    bool mesh_mapping = (*m_sim_ptr).has_mesh_mapping();

    auto& grad_p = (*m_sim_ptr).repo().get_field("gp");
    auto& pressure = (*m_sim_ptr).repo().get_field("p");
    auto& velocity = (*m_sim_ptr).repo().get_field("velocity");
    amr_wind::Field const* mesh_fac =
        mesh_mapping ? &((*m_sim_ptr)
                             .repo()
                             .get_mesh_mapping_field(amr_wind::FieldLoc::CELL))
                     : nullptr;
    amr_wind::Field const* mesh_detJ =
        mesh_mapping ? &((*m_sim_ptr)
                             .repo()
                             .get_mesh_mapping_detJ(amr_wind::FieldLoc::CELL))
                     : nullptr;
    const auto* ref_density =
        is_anelastic ? &((*m_sim_ptr).repo().get_field("reference_density"))
                     : nullptr;

    // ASA does not think we actually need to define sigma here. It
    // should not be used by calcGradPhi

    // Create sigma while accounting for mesh mapping
    // sigma = 1/(fac^2)*J * dt/rho
    int finest_level = (*m_sim_ptr).mesh().finestLevel();
    Vector<amrex::MultiFab> sigma(finest_level + 1);
    if (variable_density || mesh_mapping) {
        int ncomp = mesh_mapping ? AMREX_SPACEDIM : 1;
        for (int lev = 0; lev <= finest_level; ++lev) {
            sigma[lev].define(
                (*m_sim_ptr).mesh().boxArray(lev),
                (*m_sim_ptr).mesh().DistributionMap(lev), ncomp, 0, MFInfo());
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
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

    // Set up projection object
    std::unique_ptr<Hydro::NodalProjector> nodal_projector;

    auto bclo = nodal_projection::get_projection_bc(
        amrex::Orientation::low, pressure, (*m_sim_ptr).mesh().Geom(0));
    auto bchi = nodal_projection::get_projection_bc(
        amrex::Orientation::high, pressure, (*m_sim_ptr).mesh().Geom(0));

    // Velocity multifab is needed for proper initialization, but only the size
    // matters for the purpose of calculating gradp, the values do not matter
    Vector<MultiFab*> vel;
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel.push_back(&(velocity(lev)));
    }

    if (is_anelastic) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            amrex::Multiply(
                velocity(lev), (*ref_density)(lev), 0, 0, density[lev]->nComp(),
                0);
        }
    }

    amr_wind::MLMGOptions options("nodal_proj");

    if (variable_density || mesh_mapping) {
        nodal_projector = std::make_unique<Hydro::NodalProjector>(
            vel, GetVecOfConstPtrs(sigma),
            (*m_sim_ptr).mesh().Geom(0, finest_level), options.lpinfo());
    } else {
        amrex::Real rho_0 = 1.0;
        amrex::ParmParse pp("incflo");
        pp.query("density", rho_0);

        nodal_projector = std::make_unique<Hydro::NodalProjector>(
            vel, scaling_factor / rho_0,
            (*m_sim_ptr).mesh().Geom(0, finest_level), options.lpinfo());
    }

    // Set MLMG and NodalProjector options
    options(*nodal_projector);
    nodal_projector->setDomainBC(bclo, bchi);

    // Recalculate gradphi with fluxes
    auto gradphi = nodal_projector->calcGradPhi(pressure.vec_ptrs());

    // Transfer pressure gradient to gp field
    for (int lev = 0; lev <= finest_level; lev++) {

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(grad_p(lev), TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& tbx = mfi.tilebox();
            Array4<Real> const& gp_lev = grad_p(lev).array(mfi);
            Array4<Real const> const& gp_proj = gradphi[lev]->const_array(mfi);
            amrex::ParallelFor(
                tbx, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    gp_lev(i, j, k, n) = gp_proj(i, j, k, n);
                });
        }
    }

    // Average down is unnecessary because it is built into calcGradPhi
}

} // namespace amr_wind