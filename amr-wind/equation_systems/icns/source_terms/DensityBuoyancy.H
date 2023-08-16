//
//  DensityBuoyancy.H
//  amr-wind
//

#ifndef DensityBuoyancy_H
#define DensityBuoyancy_H

#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"

namespace amr_wind::pde::icns {

/** Buoyancy term from density variations
 *  \ingroup icns_src
 *
 *  \sa RayleighTaylor
 */
class DensityBuoyancy : public MomentumSource::Register<DensityBuoyancy>
{
public:
    static std::string identifier() { return "DensityBuoyancy"; }

    explicit DensityBuoyancy(const CFDSim& sim);

    ~DensityBuoyancy() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& vel_forces) const override;

private:
    amrex::Vector<amrex::Real> m_gravity{{0.0, 0.0, -9.81}};

    amrex::Real rho_0{1.0};

    Field& m_density;
};

} // namespace amr_wind::pde::icns

#endif /* DensityBuoyancy_H */
