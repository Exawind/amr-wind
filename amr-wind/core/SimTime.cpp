#include "amr-wind/core/SimTime.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Print.H"

namespace amr_wind {

void SimTime::parse_parameters()
{
    // Initialize deltaT to negative values
    for (amrex::Real& i : m_dt) {
        i = -1.0;
    }

    // Parse options in "time" namespace
    amrex::ParmParse pp("time");
    pp.query("stop_time", m_stop_time);
    pp.query("max_step", m_stop_time_index);
    pp.query("fixed_dt", m_fixed_dt);
    pp.query("initial_dt", m_initial_dt);
    pp.query("init_shrink", m_init_shrink);
    pp.query("cfl", m_max_cfl);
    pp.query("verbose", m_verbose);
    pp.query("regrid_interval", m_regrid_interval);
    pp.query("plot_interval", m_plt_interval);
    pp.query("checkpoint_interval", m_chkpt_interval);
    pp.query("regrid_start", m_regrid_start_index);
    pp.query("plot_start", m_plt_start_index);
    pp.query("checkpoint_start", m_chkpt_start_index);
    pp.query("use_force_cfl", m_use_force_cfl);

    if (m_fixed_dt > 0.0) {
        m_dt[0] = m_fixed_dt;
    } else {
        m_adaptive = true;
    }
}

bool SimTime::new_timestep()
{
    bool continue_sim = continue_simulation();

    if (m_is_init && (m_verbose >= 0)) {
        amrex::Print() << "\nBegin simulation: " << std::endl;
        if ((m_stop_time > 0) && (m_stop_time_index >= 0)) {
            amrex::Print() << "  Run until " << m_stop_time << " sec. or "
                           << m_stop_time_index << " timesteps" << std::endl;
        } else if (m_stop_time > 0) {
            amrex::Print() << "  Run till " << m_stop_time << " seconds "
                           << std::endl;
        } else if (m_stop_time_index >= 0) {
            amrex::Print() << "  Run for " << m_stop_time_index << " timesteps"
                           << std::endl;
        }
        if (m_adaptive) {
            amrex::Print() << "  Adaptive timestepping with max. CFL = "
                           << m_max_cfl << std::endl;
        } else {
            amrex::Print() << "  Fixed timestepping with dt = " << m_fixed_dt
                           << "; max. CFL from inputs = " << m_max_cfl
                           << std::endl;
        }
    }

    // Toggle initialization state and enter evolution phase
    m_is_init = false;
    if (continue_sim) {
        for (int i = max_time_states - 1; i > 0; --i) {
            m_dt[i] = m_dt[i - 1];
        }

        m_time_index++;
        m_cur_time = m_new_time;
        m_new_time += m_dt[0];

        // clang-format off
        if (m_verbose >= 0) {
            amrex::Print()
                << "\n==============================================================================\n";
}
        // clang-format on
    } else {
        m_cur_time = m_new_time;
    }

    return continue_sim;
}

void SimTime::set_current_cfl(
    const amrex::Real conv_cfl,
    const amrex::Real diff_cfl,
    const amrex::Real src_cfl)
{
    bool issue_cfl_warning = false;
    const amrex::Real cd_cfl = conv_cfl + diff_cfl;
    const amrex::Real cfl_unit_time =
        cd_cfl + std::sqrt(cd_cfl * cd_cfl + 4.0 * src_cfl);
    if ((m_adaptive && !m_is_init) &&
        (cfl_unit_time < std::numeric_limits<amrex::Real>::epsilon())) {
        amrex::Abort(
            "CFL is below machine epsilon and the time step is adaptive. "
            "Please use a fixed time step or fix the case setup");
    }
    amrex::Real dt_new =
        2.0 * m_max_cfl /
        amrex::max(cfl_unit_time, std::numeric_limits<amrex::Real>::epsilon());

    // Restrict timestep during initialization phase
    if (m_is_init) {
        dt_new *= m_init_shrink;
    }

    // Limit timestep growth to 10% per timestep
    if (m_dt[0] > 0.0) {
        dt_new = amrex::min<amrex::Real>(dt_new, 1.1 * m_dt[0]);
    }

    // Don't overshoot stop time
    if ((m_stop_time > 0.0) && ((m_cur_time + dt_new) > m_stop_time)) {
        dt_new = m_stop_time - m_cur_time;
    }

    if (m_adaptive) {
        m_dt[0] = dt_new;

        if (m_is_init && m_initial_dt > 0.0) {
            m_dt[0] = amrex::min(dt_new, m_initial_dt);
        }

        if (!m_is_init) {
            m_new_time = m_cur_time + m_dt[0];
        }

    } else {
        // If user has specified fixed DT then issue a warning if the timestep
        // is larger than the deltaT determined from max. CFL considerations.
        // Only issue warnings when the error is greater than 1% of the timestep
        // specified
        if ((1.0 - (dt_new / m_fixed_dt)) > 0.01) {
            issue_cfl_warning = true;
        }

        // Ensure that we use user-specified dt. Checkpoint restart might have
        // overridden this
        m_dt[0] = m_fixed_dt;
    }

    m_current_cfl = 0.5 * cfl_unit_time * m_dt[0];
    if (m_verbose >= 0) {
        if (!m_is_init) {
            amrex::Print() << "Step: " << m_time_index << " dt: " << m_dt[0]
                           << " Time: " << std::setprecision(6) << m_cur_time
                           << " to " << m_new_time << std::endl;
        } else {
            amrex::Print() << "dt: " << std::setprecision(6) << m_dt[0]
                           << std::endl;
        }
        amrex::Print() << "CFL: " << std::setprecision(6) << m_current_cfl
                       << " (conv: " << std::setprecision(6)
                       << conv_cfl * m_dt[0]
                       << " diff: " << std::setprecision(6)
                       << diff_cfl * m_dt[0] << " src: " << std::setprecision(6)
                       << std::sqrt(src_cfl) * m_dt[0] << " )" << std::endl;
    }
    if (issue_cfl_warning && !m_is_init) {
        amrex::Print() << "WARNING: fixed_dt does not satisfy CFL condition.\n"
                       << "Max. CFL: " << m_max_cfl
                       << " => dt: " << std::setprecision(6) << dt_new
                       << "; dt_inp: " << m_fixed_dt << std::endl;
    }
}

bool SimTime::continue_simulation() const
{
    constexpr double eps = 1.0e-12;
    bool stop_simulation = false;

    if (m_stop_time_index == 0) {
        return stop_simulation;
    }

    if ((m_stop_time > 0.0) && ((m_new_time + eps) >= m_stop_time)) {
        return stop_simulation;
    }

    if ((m_stop_time_index > 0) && (m_time_index >= m_stop_time_index)) {
        return stop_simulation;
    }

    return !(stop_simulation);
}

bool SimTime::do_regrid() const
{
    return (
        (m_regrid_interval > 0) &&
        ((m_time_index - m_regrid_start_index) > 0) &&
        ((m_time_index - m_regrid_start_index) % m_regrid_interval == 0));
}

bool SimTime::write_plot_file() const
{
    return (
        (m_plt_interval > 0) &&
        ((m_time_index - m_plt_start_index) % m_plt_interval == 0));
}

bool SimTime::write_checkpoint() const
{
    return (
        (m_chkpt_interval > 0) &&
        ((m_time_index - m_chkpt_start_index) % m_chkpt_interval == 0));
}

bool SimTime::write_last_plot_file() const
{
    return (
        (m_plt_interval > 0) &&
        ((m_time_index - m_plt_start_index) % m_plt_interval != 0));
}

bool SimTime::write_last_checkpoint() const
{
    return (
        (m_chkpt_interval > 0) &&
        ((m_time_index - m_chkpt_start_index) % m_chkpt_interval != 0));
}

void SimTime::set_restart_time(int tidx, amrex::Real time)
{
    m_time_index = tidx;
    m_start_time_index = tidx;

    m_new_time = time;
    m_cur_time = time;
    m_start_time = time;
}

} // namespace amr_wind
