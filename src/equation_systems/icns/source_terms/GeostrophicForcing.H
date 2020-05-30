#ifndef GEOSTROPHICFORCING_H
#define GEOSTROPHICFORCING_H

#include "icns/MomentumSource.H"
#include "SimTime.H"

namespace amr_wind {
namespace pde {
namespace icns {

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
    const Field& m_density;
    
    //! Target velocity
    amrex::Vector<amrex::Real> m_target_vel{{0.0, 0.0, 0.0}};

    //! Forcing source term (pressure gradient)
    amrex::Vector<amrex::Real> m_g_forcing{{0.0, 0.0, 0.0}};

};

} // namespace icns
} // namespace pde
} // namespace amr_wind

#endif /* GEOSTROPHICFORCING_H */
