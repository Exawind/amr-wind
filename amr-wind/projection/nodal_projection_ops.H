#ifndef NODAL_PROJECTION_OPS_H
#define NODAL_PROJECTION_OPS_H

#include <AMReX_BC_TYPES.H>
#include "amr-wind/core/field_ops.H"
#include "amr-wind/core/Physics.H"
#include "amr-wind/wind_energy/ABL.H"

using namespace amrex;

namespace amr_wind::nodal_projection {

void set_inflow_velocity(
    amr_wind::PhysicsMgr& phy_mgr,
    amr_wind::Field& vel_fld,
    int lev,
    amrex::Real time,
    MultiFab& vel_mfab,
    int nghost);

Array<amrex::LinOpBCType, AMREX_SPACEDIM> get_projection_bc(
    Orientation::Side side, amr_wind::Field& pressure, amrex::Geometry geom0);

} // namespace amr_wind::nodal_projection

#endif