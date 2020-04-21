#ifndef BOUSSINESQBUOYANCY_H
#define BOUSSINESQBUOYANCY_H

#include "AMReX_Array4.H"
#include "AMReX_Box.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"
#include "AMReX_MultiFab.H"
#include "FieldRepo.H"
#include "icns/MomentumSource.H"
#include "SimTime.H"

namespace amr_wind {

class BoussinesqBuoyancy : public pde::MomentumSource
{
public:
    BoussinesqBuoyancy(const SimTime& time_in, const FieldRepo& repo_in);

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const;

private:
    amrex::Vector<amrex::Real> m_gravity{{0.0, 0.0, -9.81}};

    //! Reference temperature (Kelvin)
    amrex::Real m_ref_theta{300.0};

    //! Thermal expansion coefficient
    amrex::Real m_beta{0.0};

    Field& m_temperature;

};

}

#endif /* BOUSSINESQBUOYANCY_H */
