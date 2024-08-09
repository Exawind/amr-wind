//
//  Godunov.H
//  amr-wind
//

#ifndef Godunov_H
#define Godunov_H

#include "amr-wind/core/FieldRepo.H"

namespace godunov {

enum class scheme { PLM, PPM, PPM_NOLIM, BDS, WENO_JS, WENOZ, MINMOD, UPWIND };

void compute_fluxes(
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
    amrex::Array4<amrex::Real const> const& fq,
    amrex::BCRec const* pbc,
    int const* iconserv,
    amrex::Real* p,
    amrex::Vector<amrex::Geometry> geom,
    amrex::Real dt,
    godunov::scheme godunov_scheme,
    bool godunov_use_forces_in_trans);

} // namespace godunov

#endif /* Godunov_H */
