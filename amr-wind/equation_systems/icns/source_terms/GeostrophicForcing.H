#ifndef GEOSTROPHICFORCING_H
#define GEOSTROPHICFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"

namespace amr_wind {
namespace pde {
namespace icns {

/** Apply source term equivalent to drive the flow to match given geostrophic wind
 *  \ingroup icns_src we_abl
 */
class GeostrophicForcing : public MomentumSource::Register<GeostrophicForcing>
{
public:
    static const std::string identifier() { return "GeostrophicForcing"; }

    GeostrophicForcing(const CFDSim&);

    virtual ~GeostrophicForcing();

    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    
    //! Target velocity
    amrex::Vector<amrex::Real> m_target_vel{{0.0, 0.0, 0.0}};

    //! Forcing source term (pressure gradient)
    amrex::Vector<amrex::Real> m_g_forcing{{0.0, 0.0, 0.0}};

};

} // namespace icns
} // namespace pde
} // namespace amr_wind

#endif /* GEOSTROPHICFORCING_H */
