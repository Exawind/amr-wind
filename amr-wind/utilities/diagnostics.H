#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include "AMReX_MultiFabUtil.H"

namespace amr_wind::diagnostics {

amrex::Real get_vel_max(
    const amrex::MultiFab& vel,
    const amrex::iMultiFab& level_mask,
    const int vdir,
    const amrex::Real factor);

amrex::Real get_vel_max(
    const amrex::MultiFab& vel,
    const amrex::iMultiFab& level_mask,
    const int vdir);

amrex::Real get_vel_min(
    const amrex::MultiFab& vel,
    const amrex::iMultiFab& level_mask,
    const int vdir);

amrex::Real get_vel_loc(
    amrex::MultiFab& vel,
    amrex::iMultiFab& level_mask,
    int vdir,
    int ldir,
    amrex::Real vel_max,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> problo,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx);

amrex::Real get_macvel_max(
    const amrex::MultiFab& vel,
    const amrex::iMultiFab& level_mask,
    const int vdir,
    const amrex::Real factor);

amrex::Real get_macvel_max(
    const amrex::MultiFab& vel,
    const amrex::iMultiFab& level_mask,
    const int vdir);

amrex::Real get_macvel_min(
    const amrex::MultiFab& vel,
    const amrex::iMultiFab& level_mask,
    const int vdir);

amrex::Real get_macvel_loc(
    amrex::MultiFab& macvel,
    amrex::iMultiFab& level_mask,
    int vdir,
    int ldir,
    amrex::Real vel_max,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> problo,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx);

amrex::Array<amrex::Real, 24> PrintMaxVelLocations(
    const amr_wind::FieldRepo& repo, const std::string& header);

amrex::Array<amrex::Real, 24> PrintMaxMACVelLocations(
    const amr_wind::FieldRepo& repo, const std::string& header);

} // namespace amr_wind::diagnostics

#endif