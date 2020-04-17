//
//  DensityBuoyancy.H
//  amr-wind
//


#ifndef DensityBuoyancy_H
#define DensityBuoyancy_H

#include "AMReX_Array4.H"
#include "AMReX_Box.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {

class DensityBuoyancy
{
public:
    DensityBuoyancy();

    void operator()(
        const amrex::Box& bx,
        const amrex::Array4<const amrex::Real>& density,
        const amrex::Array4<amrex::Real>& vel_forces) const;

private:
    amrex::Vector<amrex::Real> m_gravity{{0.0, 0.0, -9.81}};
    amrex::Real rho_0{1.0};

};

}

#endif /* DensityBuoyancy_H */
