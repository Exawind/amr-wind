#ifndef ABLFORCING_H
#define ABLFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"

namespace amr_wind::pde::icns {

/** Forcing term to drive the ABL profile to a desired state
 *  \ingroup icns_src we_abl
 *
 *  \sa ABL
 */
class ABLForcing : public MomentumSource::Register<ABLForcing>
{
public:
    static std::string identifier() { return "ABLForcing"; }

    explicit ABLForcing(const CFDSim& sim);

    ~ABLForcing() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

    inline void set_target_velocities(amrex::Real ux, amrex::Real uy)
    {
        m_target_vel[0] = ux;
        m_target_vel[1] = uy;
    }

    inline void set_mean_velocities(amrex::Real ux, amrex::Real uy)
    {
        m_mean_vel[0] = ux;
        m_mean_vel[1] = uy;

        const auto& current_time = m_time.current_time();
        const auto& dt = m_time.deltaT();

        if (!m_vel_timetable.empty()) {
            const amrex::Real current_spd = ::amr_wind::interp::linear(
                m_time_table, m_speed_table, current_time);
            const amrex::Real current_dir = ::amr_wind::interp::linear(
                m_time_table, m_direction_table, current_time);

            m_target_vel[0] = current_spd * std::cos(current_dir);
            m_target_vel[1] = current_spd * std::sin(current_dir);
        }

        m_abl_forcing[0] = (m_target_vel[0] - m_mean_vel[0]) / dt;
        m_abl_forcing[1] = (m_target_vel[1] - m_mean_vel[1]) / dt;
    }

    amrex::RealArray abl_forcing() const { return m_abl_forcing; }

    amrex::Real forcing_height() const { return m_forcing_height; }

private:
    const SimTime& m_time;

    //!  ABL forcing terms
    amrex::RealArray m_abl_forcing{{0.0, 0.0, 0.0}};

    //! File name for velocity forcing time table
    std::string m_vel_timetable;

    //! Velocity forcing time table
    amrex::Vector<amrex::Real> m_time_table;

    //! Velocity forcing speed table
    amrex::Vector<amrex::Real> m_speed_table;

    //! Velocity forcing direction table
    amrex::Vector<amrex::Real> m_direction_table;

    //! Target velocity
    amrex::Vector<amrex::Real> m_target_vel{{0.0, 0.0, 0.0}};

    //! Current mean vel
    amrex::RealArray m_mean_vel{{0.0, 0.0, 0.0}};

    //! Height at which the velocities are forcing
    amrex::Real m_forcing_height;
};

} // namespace amr_wind::pde::icns

#endif /* ABLFORCING_H */
