#ifndef BOUSSINESQBUOYANCY_H
#define BOUSSINESQBUOYANCY_H

#include "AMReX_Array4.H"
#include "AMReX_Box.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {

class BoussinesqBuoyancy
{
public:
    BoussinesqBuoyancy();

    void operator()(
        const amrex::Box& bx,
        const amrex::Array4<const amrex::Real>& scalars,
        const amrex::Array4<amrex::Real>& vel_forces) const;

private:
    amrex::Vector<amrex::Real> m_gravity{{0.0, 0.0, -9.81}};

    //! Reference temperature (Kelvin)
    amrex::Real m_ref_theta{300.0};

    //! Thermal expansion coefficient
    amrex::Real m_beta{0.0};
};

}

#endif /* BOUSSINESQBUOYANCY_H */
