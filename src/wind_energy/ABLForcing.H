#ifndef ABLFORCING_H
#define ABLFORCING_H

#include "icns/MomentumSource.H"
#include "SimTime.H"

class PlaneAveraging;

namespace amr_wind {

class ABLForcingOld
{
public:
    ABLForcingOld(const SimTime& time);

    void operator()(
        const amrex::Box& vbx,
        const amrex::Array4<amrex::Real>& vel_forces) const;

    inline void set_mean_velocities(amrex::Real ux, amrex::Real uy)
    {
        m_mean_vel[0] = ux;
        m_mean_vel[1] = uy;
    }

    amrex::Real forcing_height() const { return m_forcing_height; }

    const amrex::RealArray& target_velocity() const { return m_target_vel; }

private:
    const SimTime& m_time;

    //! Target velocity
    amrex::RealArray m_target_vel{{0.0, 0.0, 0.0}};

    //! Current mean vel
    amrex::RealArray m_mean_vel{{0.0, 0.0, 0.0}};

    //! Height at which the velocities are forcing
    amrex::Real m_forcing_height;
};

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
    amrex::RealArray m_target_vel{{0.0, 0.0, 0.0}};

    //! Current mean vel
    amrex::RealArray m_mean_vel{{0.0, 0.0, 0.0}};

    //! Height at which the velocities are forcing
    amrex::Real m_forcing_height;
};

} // namespace icns
} // namespace pde
} // namespace amr_wind

#endif /* ABLFORCING_H */
