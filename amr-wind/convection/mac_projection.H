//
//  mac_projection.H
//  amr-wind
//

#ifndef mac_projection_h
#define mac_projection_h

#include "amr-wind/core/FieldRepo.H"

namespace mac {

amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> get_projection_bc(
    amrex::Orientation::Side side,
    amrex::GpuArray<BC, AMREX_SPACEDIM * 2> bctype,
    amrex::Vector<amrex::Geometry> geom) noexcept;

} // namespace mac

#endif /* mac_projection_h */
