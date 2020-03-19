#ifndef ABLFORCING_H
#define ABLFORCING_H

#include "AMReX_Array.H"
#include "AMReX_Array4.H"
#include "AMReX_Box.H"
#include "AMReX_Geometry.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"

#include "SimTime.H"

namespace amr_wind {

class ABLForcing
{
public:
    ABLForcing(const SimTime& time);

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
    amrex::RealArray m_target_vel{0.0, 0.0, 0.0};

    //! Current mean vel
    amrex::RealArray m_mean_vel{0.0, 0.0, 0.0};

    //! Height at which the velocities are forcing
    amrex::Real m_forcing_height;
};

}

#endif /* ABLFORCING_H */
