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

#include "icns/MomentumSource.H"
#include "SimTime.H"

namespace amr_wind {

class DensityBuoyancy : public pde::MomentumSource
{
    
public:
    DensityBuoyancy(const SimTime& time_in, const FieldRepo& repo_in);

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const;

private:

    amrex::Vector<amrex::Real> m_gravity{{0.0, 0.0, -9.81}};
    amrex::Real rho_0{1.0};

    Field& m_density;

};

} // namespace amr_wind

#endif /* DensityBuoyancy_hpp */
