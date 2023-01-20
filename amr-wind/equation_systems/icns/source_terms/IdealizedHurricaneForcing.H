#ifndef IDEALIZED_HURRICANE_FORCING_H
#define IDEALIZED_HURRICANE_FORCING_H

#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"

namespace amr_wind::pde::icns {

/** Apply mesoscale tendency source term to drive the flow
 * and match the hurricane cyclic/rotational forcing at large
 * distances from the eye-wall. This implements the "Simple"
 * method of George Bryant.
 * Details can be found in:
 * "A simple Method for Simulating Wind Profiles in the Boundary Layer
 * of Tropical Cyclones. Bryant et al. Boundary-Layer Meteorol (2017)
 * DOI 10.1007/s10546-016-0207-0
 */

class IdealizedHurricaneForcing
    : public MomentumSource::Register<IdealizedHurricaneForcing>
{
public:
    static std::string identifier() { return "IdealizedHurricaneForcing"; }

    explicit IdealizedHurricaneForcing(const CFDSim& /*unused*/);

    ~IdealizedHurricaneForcing() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    //! Velocity field
    const Field& m_velocity;

    //! Gradient velocity
    amrex::Real m_Ug{25.0}; // Ug=25 m/s is the default value

    //! Eye-wall radial distance
    amrex::Real m_R{40000.0}; // 40km is taken as the default value
};

} // namespace amr_wind::pde::icns

#endif /* IDEALIZED_HURRICANE_FORCING_H */