#ifndef DamBreakFIELDINIT_H
#define DamBreakFIELDINIT_H

#include "AMReX_Array.H"
#include "AMReX_Array4.H"
#include "AMReX_Box.H"
#include "AMReX_Geometry.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"
#include "AMReX_Gpu.H"

namespace amr_wind {

/** Initialize mutiphase DamBreak fields
 */

class DamBreakFieldInit
{
    static_assert(AMREX_SPACEDIM == 3, "DamBreak requires 3 dimensional mesh");

public:
    DamBreakFieldInit();

    void operator()(
        const amrex::Box& vbx,
        const amrex::Geometry& geom,
        const amrex::Array4<amrex::Real>& velocity,
        const amrex::Array4<amrex::Real>& levelset) const;

private:
    //! Initial DamBreak location
    amrex::Vector<amrex::Real> m_loc{{0.0, 0.0, 0.0}};

    //! dam width value
    amrex::Real m_width{0.1};

    //! dam height value
    amrex::Real m_height{0.25};
};

} // namespace amr_wind

#endif /* DamBreakFIELDINIT_H */