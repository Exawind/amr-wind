//
//  mac_projection.H
//  amr-wind
//

#ifndef mac_projection_h
#define mac_projection_h

#include "FieldRepo.H"

namespace mac {

amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> get_projection_bc(
    amrex::Orientation::Side side,
    amrex::GpuArray<BC, AMREX_SPACEDIM * 2> bctype,
    amrex::Vector<amrex::Geometry> geom) noexcept;

void apply_MAC_projection(
    amr_wind::FieldRepo& repo,
    amr_wind::FieldState fstate,
    int max_coarsening_level,
    amrex::Real mac_mg_rtol,
    amrex::Real mac_mg_atol);

} // namespace mac

#endif /* mac_projection_h */
