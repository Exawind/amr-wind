#ifndef SplitAdvection_H
#define SplitAdvection_H

#include "amr-wind/core/FieldRepo.H"

namespace amr_wind {
namespace multiphase {

void split_advection(
    int lev,
    amrex::Box const& bx,
    int isweep,
    amrex::Array4<amrex::Real> const& volfrac,
    amrex::Array4<amrex::Real const> const& umac,
    amrex::Array4<amrex::Real const> const& vmac,
    amrex::Array4<amrex::Real const> const& wmac,
    amrex::BCRec const* pbc,
    amrex::Real* p,
    amrex::Vector<amrex::Geometry> geom,
    amrex::Real dt,
    bool is_lagrangian);

void sweep(
    const int dir,
    amrex::Box const& bx,
    amrex::Real dtdx,
    amrex::Array4<amrex::Real const> const& vel_mac,
    amrex::Array4<amrex::Real> const& volfrac,
    amrex::Array4<amrex::Real> const& fluxL,
    amrex::Array4<amrex::Real> const& fluxC,
    amrex::Array4<amrex::Real> const& fluxR,
    amrex::BCRec const* pbc,
    const int dimLow,
    const int dimHigh,
    bool is_lagrangian);

} // namespace multiphase
} // namespace amr_wind
#endif /* SplitAdvection_H */