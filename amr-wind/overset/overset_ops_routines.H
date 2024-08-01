#ifndef OVERSET_OPS_ROUTINES_H_
#define OVERSET_OPS_ROUTINES_H_

#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFab.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "amr-wind/overset/overset_ops_K.H"

namespace amr_wind::overset_ops {

// Populate approximate signed distance function using vof field
void populate_psi(
    amrex::MultiFab& mf_psi,
    const amrex::MultiFab& mf_vof,
    const amrex::Real i_th,
    const amrex::Real asdf_tiny);

// Modify a vof field to not have values that barely differ from 0 or 1
void process_vof(amrex::MultiFab& mf_vof, const amrex::Real vof_tol);

// Combine overset target vof field with current non-overset vof field
void harmonize_vof(
    amrex::MultiFab& mf_vof_target,
    const amrex::MultiFab& mf_vof_original,
    const amrex::iMultiFab& mf_iblank);

// Populate normal vector with special treatment of overset boundary
void populate_normal_vector(
    amrex::MultiFab& mf_normvec,
    const amrex::MultiFab& mf_vof,
    const amrex::iMultiFab& mf_iblank);

// Calculate fluxes for reinitialization over entire domain without concern for
// overset bdy
void populate_sharpen_fluxes(
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    const amrex::MultiFab& mf_vof,
    const amrex::MultiFab& mf_target_vof,
    const amrex::MultiFab& mf_norm,
    const amrex::MultiFab& mf_velocity,
    const amrex::MultiFab& mf_gp,
    const amrex::MultiFab& mf_density,
    const amrex::Real Gamma,
    const amrex::Real margin,
    const amrex::Real rho1,
    const amrex::Real rho2);

// Process reinitialization fluxes - zero non-internal to overset region;
// also calculate pressure source / sink term as a function of fluxes
void process_fluxes_calc_src(
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    amrex::MultiFab& mf_psource,
    const amrex::MultiFab& mf_vof,
    const amrex::iMultiFab& mf_iblank);

amrex::Real calculate_pseudo_velocity_scale(
    const amrex::iMultiFab& mf_iblank,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
    const amrex::Real pvmax);

// Calculate a type of CFL by measuring how much % VOF is being removed per cell
amrex::Real calculate_pseudo_dt_flux(
    const amrex::MultiFab& mf_fx,
    const amrex::MultiFab& mf_fy,
    const amrex::MultiFab& mf_fz,
    const amrex::MultiFab& mf_vof,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
    const amrex::Real tol);

// Apply reinitialization fluxes to modify fields
void apply_fluxes(
    const amrex::MultiFab& mf_fx,
    const amrex::MultiFab& mf_fy,
    const amrex::MultiFab& mf_fz,
    const amrex::MultiFab& mf_psource,
    amrex::MultiFab& mf_vof,
    amrex::MultiFab& mf_dens,
    amrex::MultiFab& mf_vel,
    amrex::MultiFab& mf_gp,
    amrex::MultiFab& mf_pressure,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
    const amrex::Real ptfac,
    const amrex::Real vof_tol);

// Get the size of the smallest VOF flux to quantify convergence
amrex::Real measure_convergence(
    amrex::MultiFab& mf_fx, amrex::MultiFab& mf_fy, amrex::MultiFab& mf_fz);

// Set levelset field to another quantity to view in plotfile for debugging
void equate_field(amrex::MultiFab& mf_dest, const amrex::MultiFab& mf_src);

// Replace pressure gradient with hydrostatic field in overset regions
void replace_gradp_hydrostatic(
    amrex::MultiFab& mf_gp,
    const amrex::MultiFab& mf_density,
    const amrex::MultiFab& mf_refdens,
    const amrex::iMultiFab& mf_iblank,
    const amrex::Real grav_z,
    const bool is_pptb);

// Swap pressure gradient values in overset region
void replace_gradp(
    amrex::MultiFab& mf_gp,
    const amrex::MultiFab& mf_gp0,
    const amrex::iMultiFab& mf_iblank);

// Apply pressure gradient to velocity field
void apply_pressure_gradient(
    amrex::MultiFab& mf_vel,
    const amrex::MultiFab& mf_density,
    const amrex::MultiFab& mf_gp,
    const amrex::Real scaling_factor);

} // namespace amr_wind::overset_ops

#endif