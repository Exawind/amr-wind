#ifndef SplitLagrangianAdvection_H
#define SplitLagrangianAdvection_H

#include "amr-wind/core/FieldRepo.H"

namespace amr_wind {
namespace multiphase {

void split_lagrangian_advection(
    int lev,
    amrex::Box const& bx,
    amrex::Array4<amrex::Real const> const& q,
    amrex::Array4<amrex::Real const> const& umac,
    amrex::Array4<amrex::Real const> const& vmac,
    amrex::Array4<amrex::Real const> const& wmac,
    amrex::BCRec const* pbc,
    amrex::Vector<amrex::Geometry> geom,
    amrex::Real dt);

} // namespace multiphase
} // namespace amr_wind
#endif /* SplitLagrangianAdvection_H */