#ifndef GEOSTROPHICFORCING_H
#define GEOSTROPHICFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"

namespace amr_wind::pde::icns {

/** Apply source term equivalent to drive the flow to match given geostrophic
 * wind.
 *
 * Geostrophic wind forcing source term is the cross product of the Coriolis
 * factor with the Geostrophic wind velocity field vector. Rotating the Coriolis
 * factor into the local frame of reference from the earths frame of reference
 * results in two components Coriolis_factor(cosphi, sinphi) in the y and z
 * directions respectively. Taking the cross product of Coriolis_factor(cosphi,
 * sinphi) with the Geostrophic wind velocity field vector U_g(x,y,z) results in
 * three source term components src(x,y,z). The vertical component of the
 * Geostrophic wind source term is typically neglected as it is unphysical and
 * is also, negligble compared to the other vertical component source terms.
 * There is an option to include the vertical Geostrophic source term by turning
 * on three_ComponentForcing = true. three_ComponentForcing is by default set to
 * false. This removes the contributions of the Coriolis_factor(cosphi)
 * component from the Geostrophic wind source term components.
 *
 * \ingroup icns_src we_abl
 */
class GeostrophicForcing : public MomentumSource::Register<GeostrophicForcing>
{
public:
    static std::string identifier() { return "GeostrophicForcing"; }

    explicit GeostrophicForcing(const CFDSim& /*unused*/);

    ~GeostrophicForcing() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    //! Latitude
    amrex::Real latitude = {0.0};
    amrex::Real sinphi{0.0};
    amrex::Real cosphi{0.0};

    //! Default three component forcing
    bool m_three_dimensional_forcing{false};

    //! Target velocity
    amrex::Vector<amrex::Real> m_target_vel{{0.0, 0.0, 0.0}};

    //! Forcing source term (pressure gradient)
    amrex::Vector<amrex::Real> m_g_forcing{{0.0, 0.0, 0.0}};
};

} // namespace amr_wind::pde::icns

#endif /* GEOSTROPHICFORCING_H */
