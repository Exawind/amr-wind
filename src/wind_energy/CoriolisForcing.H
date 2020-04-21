#ifndef CORIOLISFORCING_H
#define CORIOLISFORCING_H

#include "AMReX_Array4.H"
#include "AMReX_Box.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"
#include "AMReX_MultiFab.H"

#include "FieldRepo.H"
#include "icns/MomentumSource.H"
#include "SimTime.H"

namespace amr_wind {

class CoriolisForcing : public pde::MomentumSource
{
public:
    CoriolisForcing(const SimTime& time_in, const FieldRepo& repo_in);

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
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

    Field& m_velocity;
    
};

}

#endif /* CORIOLISFORCING_H */
