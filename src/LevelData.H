#ifndef LEVELDATA_H
#define LEVELDATA_H

#include "AMReX_AmrCore.H"
#include "AMReX_MultiFab.H"

struct LevelData
{
    LevelData() = default;
    LevelData(
        amrex::BoxArray const& ba,
        amrex::DistributionMapping const& dm,
        amrex::FabFactory<amrex::FArrayBox> const& fact,
        int ntrac,
        int ng_state,
        bool use_godunov,
        bool implicit_diffusion,
        bool advect_tracer);
    // cell-centered multifabs
    amrex::MultiFab velocity;
    amrex::MultiFab velocity_o;
    amrex::MultiFab density;
    amrex::MultiFab density_o;
    amrex::MultiFab tracer;
    amrex::MultiFab tracer_o;
    amrex::MultiFab gp; // pressure gradient
    // nodal multifabs
    amrex::MultiFab p;
    //
    amrex::MultiFab conv_velocity;
    amrex::MultiFab conv_velocity_o;
    amrex::MultiFab conv_density;
    amrex::MultiFab conv_density_o;
    amrex::MultiFab conv_tracer;
    amrex::MultiFab conv_tracer_o;
    //
    amrex::MultiFab divtau;
    amrex::MultiFab divtau_o;
    amrex::MultiFab laps;
    amrex::MultiFab laps_o;
};

#endif /* LEVELDATA_H */
