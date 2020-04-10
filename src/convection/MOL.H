//
//  MOL.H
//  amr-wind
//
//

#ifndef MOL_H
#define MOL_H

#include "FieldRepo.H"

namespace mol {

void compute_convective_term(
    amr_wind::FieldRepo& repo,
    const amr_wind::FieldState fstate,
    bool constant_density,
    bool advect_tracer);

void compute_convective_fluxes(
    int lev,
    amrex::Box const& bx,
    int ncomp,
    amrex::Array4<amrex::Real> const& fx,
    amrex::Array4<amrex::Real> const& fy,
    amrex::Array4<amrex::Real> const& fz,
    amrex::Array4<amrex::Real const> const& q,
    amrex::Array4<amrex::Real const> const& umac,
    amrex::Array4<amrex::Real const> const& vmac,
    amrex::Array4<amrex::Real const> const& wmac,
    amrex::BCRec const* h_bcrec,
    amrex::BCRec const* d_bcrec,
    amrex::Vector<amrex::Geometry> geom);

void compute_convective_rate(
    amrex::Box const& bx,
    int ncomp,
    amrex::Array4<amrex::Real> const& dUdt,
    amrex::Array4<amrex::Real const> const& fx,
    amrex::Array4<amrex::Real const> const& fy,
    amrex::Array4<amrex::Real const> const& fz,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxi);

} // namespace mol

#endif /* MOL_H */
