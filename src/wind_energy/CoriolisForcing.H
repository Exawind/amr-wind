#ifndef CORIOLISFORCING_H
#define CORIOLISFORCING_H

#include "AMReX_Array4.H"
#include "AMReX_Box.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"

namespace amr_wind {

class CoriolisForcingOld
{
public:
    CoriolisForcingOld();

    void operator()(
        const amrex::Box& bx,
        const amrex::Array4<amrex::Real const>& velocity,
        const amrex::Array4<amrex::Real>& vel_forces) const;

private:
    ///@{
    /** Orientation of the grid w.r.t. to planetary coordinate system
     */
    amrex::Vector<amrex::Real> m_east{{1.0, 0.0, 0.0}};
    amrex::Vector<amrex::Real> m_north{{0.0, 1.0, 0.0}};
    amrex::Vector<amrex::Real> m_up{{1.0, 0.0, 1.0}};
    ///@}

    //! Latitude where the Coriolis forcing is computed (+ve northern hemisphere)
    amrex::Real m_latitude{0.0};

    amrex::Real m_sinphi{0.0};
    amrex::Real m_cosphi{0.0};

    //! `2.0 * \Omega`
    amrex::Real m_coriolis_factor{0.0};
};

}

#endif /* CORIOLISFORCING_H */
