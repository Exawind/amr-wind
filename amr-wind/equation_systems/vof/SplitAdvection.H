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
    amrex::Real dt);

} // namespace multiphase
} // namespace amr_wind
#endif /* SplitAdvection_H */