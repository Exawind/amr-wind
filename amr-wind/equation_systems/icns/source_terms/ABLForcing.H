#ifndef ABLFORCING_H
#define ABLFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"

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

    void mean_velocity_init(const VelPlaneAveraging& /*vavg*/);

    void mean_velocity_update(const VelPlaneAveraging& /*vavg*/);

    inline void set_target_velocities(amrex::Real ux, amrex::Real uy)
    {
        m_target_vel[0] = ux;
        m_target_vel[1] = uy;
    }

    inline void set_boundary_layer_height(amrex::Real zi) { m_zi = zi; }

    inline void set_mean_velocities(amrex::Real ux, amrex::Real uy)
    {
        m_mean_vel[0] = ux;
        m_mean_vel[1] = uy;

        const auto& current_time = m_time.current_time();
        const auto& new_time = m_time.new_time();
        const auto& nph_time = 0.5 * (current_time + new_time);
        const auto& dt = m_time.delta_t();
        const auto& t_step = m_time.time_index();

        if (!m_vel_timetable.empty()) {
            // Forces should be applied at n+1/2. Because ABL forcing is a
            // delta, the difference between the target velocity (at n+1) and
            // the current velocity (at n) puts the force term at n+1/2
            const amrex::Real new_spd = ::amr_wind::interp::linear(
                m_time_table, m_speed_table, new_time);
            const amrex::Real new_dir = ::amr_wind::interp::linear_angle(
                m_time_table, m_direction_table, new_time, 2.0 * M_PI);

            m_target_vel[0] = new_spd * std::cos(new_dir);
            m_target_vel[1] = new_spd * std::sin(new_dir);
        }

        m_abl_forcing[0] = (m_target_vel[0] - m_mean_vel[0]) / dt;
        m_abl_forcing[1] = (m_target_vel[1] - m_mean_vel[1]) / dt;

        if (m_write_force_timetable &&
            amrex::ParallelDescriptor::IOProcessor() &&
            (t_step % m_force_outfreq == 0) &&
            (current_time >= m_force_outstart)) {
            std::ofstream outfile;
            // Forces are recorded at n+1/2
            outfile.open(m_force_timetable, std::ios::out | std::ios::app);
            outfile << std::setprecision(17) << nph_time << "\t"
                    << m_abl_forcing[0] << "\t" << m_abl_forcing[1] << "\t"
                    << 0.0 << "\t" << m_abl_forcing[1] / m_coriolis_factor
                    << "\t" << -m_abl_forcing[0] / m_coriolis_factor << "\t"
                    << m_zi << std::endl;
        }
    }

    amrex::RealArray abl_forcing() const { return m_abl_forcing; }

    amrex::Real forcing_height() const { return m_forcing_height; }

private:
    const SimTime& m_time;
    const amrex::AmrCore& m_mesh;

    amrex::Gpu::DeviceVector<amrex::Real> m_vel_ht;
    amrex::Gpu::DeviceVector<amrex::Real> m_vel_vals;

    //! Axis over which averages are computed.
    int m_axis{2};

    //! Activated when water is present in domain
    bool m_use_phase_ramp{false};

    //! Activate free-atmosphere damping
    bool m_fa_damping{false};

    //! Detect free-atmosphere base height (set to zi from ABLStats).
    bool m_fa_detect_height{false};

    //! Free-atmosphere base height.
    amrex::Real m_fa_height{750.0};

    //! Start time for free-atmosphere damping.
    amrex::Real m_fa_time_start{0.0};

    //! End time for free-atmosphere damping.
    amrex::Real m_fa_time_end{1.0E9};

    //! Time scale for free-atmosphere damping.
    amrex::Real m_fa_tau{100.0};

    //! Coriolis factor used for free-atmosphere damping.
    amrex::Real m_coriolis_factor;

    //! Number of cells in band to prevent forcing near liquid
    int m_n_band{2};

    //! ABL forcing terms
    amrex::RealArray m_abl_forcing{0.0};

    //! File name for velocity forcing time table
    std::string m_vel_timetable;

    //! Bool for writing forcing time table
    bool m_write_force_timetable{false};

    //! File name for forcing time table output
    std::string m_force_timetable;

    //! Output frequency for forces
    int m_force_outfreq{1};

    //! Output start time for force
    amrex::Real m_force_outstart{0.0};

    //! Velocity forcing time table
    amrex::Vector<amrex::Real> m_time_table;

    //! Velocity forcing speed table
    amrex::Vector<amrex::Real> m_speed_table;

    //! Velocity forcing direction table
    amrex::Vector<amrex::Real> m_direction_table;

    //! Target velocity
    amrex::Vector<amrex::Real> m_target_vel{0.0, 0.0, 0.0};

    //! Current mean vel
    amrex::RealArray m_mean_vel{0.0};

    //! Height at which the velocities are forcing
    amrex::Real m_forcing_height;

    //! Height from water interface where force is off
    amrex::Real m_forcing_mphase0;
    //! Height from water interface over which force is ramped
    amrex::Real m_forcing_mphase1;

    //! Local storage of interface location
    amrex::Real m_water_level;

    //! boundary layer height.
    amrex::Real m_zi;

    //! VOF field, to avoid forcing on liquid above force-off height
    const Field* m_vof;
};

} // namespace amr_wind::pde::icns

#endif /* ABLFORCING_H */
