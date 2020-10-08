#ifndef VortexPatchFIELDINIT_H
#define VortexPatchFIELDINIT_H

#include "AMReX_Array.H"
#include "AMReX_Array4.H"
#include "AMReX_Box.H"
#include "AMReX_Geometry.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"
#include "AMReX_Gpu.H"

namespace amr_wind {

/** Initialize VortexPatch fields
 */

class VortexPatchFieldInit
{
    static_assert(AMREX_SPACEDIM == 3, "VortexPatch requires 3 dimensional mesh");

public: 
    VortexPatchFieldInit();

    void operator()(
        const amrex::Box& vbx,
        const amrex::Geometry& geom,
        const amrex::Array4<amrex::Real>& velocity,
        const amrex::Array4<amrex::Real>& density,
        const amrex::Array4<amrex::Real>& levelset) const;

private:
    //! Initial VortexPatch location
    amrex::Vector<amrex::Real> m_loc{{0.35, 0.35, 0.35}};

    //! Initial density field
    amrex::Real m_rho{1.0};

    //! tracer bubble radius value
    amrex::Real m_radius{0.15};

};

}

#endif /* VortexPatchFIELDINIT_H */