#ifndef SplitAdvection_H
#define SplitAdvection_H

#include "amr-wind/core/FieldRepo.H"

namespace amr_wind {
namespace multiphase {

void split_advection_step(
    int isweep,
    int iorder,
    int nlevels,
    Field& dof_field,
    amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> const& fluxes,
    ScratchField& fluxC,
    amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> const&
        advalphas,
    Field const& u_mac,
    Field const& v_mac,
    Field const& w_mac,
    amrex::BCRec const* pbc,
    amrex::Vector<amrex::Geometry> geom,
    amrex::Real dt,
    bool use_lagrangian,
    bool rm_debris);

void split_compute_fluxes(
    const int lev,
    amrex::Box const& bx,
    const int isweep,
    amrex::Array4<amrex::Real const> const& volfrac,
    amrex::Array4<amrex::Real> const& fluxC,
    amrex::Array4<amrex::Real const> const& umac,
    amrex::Array4<amrex::Real const> const& vmac,
    amrex::Array4<amrex::Real const> const& wmac,
    amrex::Array4<amrex::Real> const& aax,
    amrex::Array4<amrex::Real> const& aay,
    amrex::Array4<amrex::Real> const& aaz,
    amrex::Array4<amrex::Real> const& fx,
    amrex::Array4<amrex::Real> const& fy,
    amrex::Array4<amrex::Real> const& fz,
    amrex::BCRec const* pbc,
    amrex::Real* p,
    amrex::Vector<amrex::Geometry> geom,
    const amrex::Real dt,
    const bool is_lagrangian);

void split_compute_sum(
    int lev,
    amrex::Box const& bx,
    int isweep,
    amrex::Array4<amrex::Real> const& volfrac,
    amrex::Array4<amrex::Real const> const& fluxC,
    amrex::Array4<amrex::Real const> const& umac,
    amrex::Array4<amrex::Real const> const& vmac,
    amrex::Array4<amrex::Real const> const& wmac,
    amrex::Array4<amrex::Real const> const& fx,
    amrex::Array4<amrex::Real const> const& fy,
    amrex::Array4<amrex::Real const> const& fz,
    amrex::Vector<amrex::Geometry> geom,
    amrex::Real dt,
    bool is_lagrangian);

void cmask_loop(
    amrex::Box const& bx,
    amrex::Array4<amrex::Real> const& volfrac,
    amrex::Array4<amrex::Real> const& fluxC,
    bool is_lagrangian);

void debris_loop(
    amrex::Box const& bx, amrex::Array4<amrex::Real> const& volfrac);

void sweep_fluxes(
    const int dir,
    amrex::Box const& bx,
    const amrex::Real dtdx,
    amrex::Array4<amrex::Real const> const& vel_mac,
    amrex::Array4<amrex::Real const> const& volfrac,
    amrex::Array4<amrex::Real> const& fluxL,
    amrex::Array4<amrex::Real> const& fluxC,
    amrex::Array4<amrex::Real> const& fluxR,
    amrex::Array4<amrex::Real> const& vofL,
    amrex::Array4<amrex::Real> const& vofR,
    const bool is_lagrangian);

void sweep_balance(
    const int dir,
    const amrex::Real dt,
    const amrex::Real dxi,
    amrex::Box const& bx,
    amrex::Array4<amrex::Real const> const& vel_mac,
    amrex::Array4<amrex::Real> const& volfrac,
    amrex::Array4<amrex::Real const> const& fluxF,
    amrex::Array4<amrex::Real const> const& fluxC,
    const bool is_lagrangian);

} // namespace multiphase
} // namespace amr_wind
#endif /* SplitAdvection_H */
