#ifndef ZalesakDiskFIELDINIT_H
#define ZalesakDiskFIELDINIT_H

#include "AMReX_Array.H"
#include "AMReX_Array4.H"
#include "AMReX_Box.H"
#include "AMReX_Geometry.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"
#include "AMReX_Gpu.H"

namespace amr_wind {

/** Initialize mutiphase ZalesakDisk fields
 */

class ZalesakDiskFieldInit
{
    static_assert(
        AMREX_SPACEDIM == 3, "ZalesakDisk requires 3 dimensional mesh");

public:
    ZalesakDiskFieldInit();

    void operator()(
        const amrex::Box& vbx,
        const amrex::Geometry& geom,
        const amrex::Array4<amrex::Real>& velocity,
        const amrex::Array4<amrex::Real>& levelset) const;

private:
    //! Initial ZalesakDisk location
    amrex::Vector<amrex::Real> m_loc{{0.5, 0.72, 0.24}};

    //! sphere radius value
    amrex::Real m_radius{0.16};

    //! slot width value
    amrex::Real m_width{0.04};

    //! slot depth
    amrex::Real m_depth{0.2};

    //! prescribed velocity period
    amrex::Real m_TT{6.0};
};

} // namespace amr_wind

#endif /* ZalesakDiskFIELDINIT_H */
