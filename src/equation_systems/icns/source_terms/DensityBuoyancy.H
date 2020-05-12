//
//  DensityBuoyancy.H
//  amr-wind
//

#ifndef DensityBuoyancy_H
#define DensityBuoyancy_H

#include "FieldRepo.H"
#include "icns/MomentumSource.H"

namespace amr_wind {
namespace pde {
namespace icns {

class DensityBuoyancy : public MomentumSource::Register<DensityBuoyancy>
{
public:
    static const std::string identifier() { return "DensityBuoyancy"; }

    DensityBuoyancy(const CFDSim& sim);

    virtual ~DensityBuoyancy();

    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    amrex::Vector<amrex::Real> m_gravity{{0.0, 0.0, -9.81}};

    amrex::Real rho_0{1.0};

    Field& m_density;
};

} // namespace icns
} // namespace pde
} // namespace amr_wind

#endif /* DensityBuoyancy_H */
