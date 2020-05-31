#ifndef ABLFORCING_H
#define ABLFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"

class PlaneAveraging;

namespace amr_wind {
namespace pde {
namespace icns {

class ABLForcing : public MomentumSource::Register<ABLForcing>
{
public:
    static const std::string identifier() { return "ABLForcing"; }

    ABLForcing(const CFDSim& sim);

    virtual ~ABLForcing();

    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

    inline void set_mean_velocities(amrex::Real ux, amrex::Real uy)
    {
        m_mean_vel[0] = ux;
        m_mean_vel[1] = uy;
    }

    amrex::Real forcing_height() const { return m_forcing_height; }

private:
    const SimTime& m_time;

    //! Target velocity
    amrex::Vector<amrex::Real> m_target_vel{{0.0, 0.0, 0.0}};

    //! Current mean vel
    amrex::RealArray m_mean_vel{{0.0, 0.0, 0.0}};

    //! Height at which the velocities are forcing
    amrex::Real m_forcing_height;
};

} // namespace icns
} // namespace pde
} // namespace amr_wind

#endif /* ABLFORCING_H */
