#include "amr-wind/overset/OversetOps.H"
#include "amr-wind/core/field_ops.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFabUtil.H"
#include "amr-wind/core/MLMGOptions.H"
#include "amr-wind/projection/nodal_projection_ops.H"
#include <hydro_NodalProjector.H>

namespace amr_wind {

void OversetOps::initialize(CFDSim& sim) { m_sim_ptr = &sim; }

void OversetOps::pre_advance_work()
{
    // Update pressure gradient using updated overset pressure field
    update_gradp();
}

void OversetOps::update_gradp()
{
    BL_PROFILE("amr-wind::OversetOps::update_gradp");

    // Get relevant fields
    auto& grad_p = (*m_sim_ptr).repo().get_field("gp");
    auto& pressure = (*m_sim_ptr).repo().get_field("p");
    auto& velocity = (*m_sim_ptr).repo().get_field("velocity");

    // Set up projection object
    std::unique_ptr<Hydro::NodalProjector> nodal_projector;
    // Boundary conditions from field, periodicity
    auto bclo = nodal_projection::get_projection_bc(
        amrex::Orientation::low, pressure,
        (*m_sim_ptr).mesh().Geom(0).isPeriodic());
    auto bchi = nodal_projection::get_projection_bc(
        amrex::Orientation::high, pressure,
        (*m_sim_ptr).mesh().Geom(0).isPeriodic());
    // Velocity multifab is needed for proper initialization, but only the size
    // matters for the purpose of calculating gradp, the values do not matter
    int finest_level = (*m_sim_ptr).mesh().finestLevel();
    Vector<MultiFab*> vel;
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel.push_back(&(velocity(lev)));
    }
    amr_wind::MLMGOptions options("nodal_proj");
    // Create nodal projector with unity scaling factor for simplicity
    nodal_projector = std::make_unique<Hydro::NodalProjector>(
        vel, 1.0, (*m_sim_ptr).mesh().Geom(0, finest_level), options.lpinfo());
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

    // Averaging down here would be unnecessary; it is built into calcGradPhi
}

} // namespace amr_wind