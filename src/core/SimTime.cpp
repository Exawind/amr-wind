#include "amr-wind/core/SimTime.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Print.H"

namespace amr_wind {

void SimTime::parse_parameters()
{
    // Initialize deltaT to negative values
    for (int i = 0; i < max_time_states; ++i) {
        m_dt[i] = -1.0;
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

    if (m_is_init) {
        amrex::Print() << "\nBegin simulation: ";
        if (m_stop_time > 0)
            amrex::Print() << "run till " << m_stop_time << " seconds" << std::endl;
        else if (m_stop_time_index >= 0)
            amrex::Print() << "run for " << m_stop_time_index << " timesteps" << std::endl;
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
        if (m_verbose >= 0)
            amrex::Print()
                << "\n==============================================================================\n"
                << "Step: " << m_time_index << " Time: " << m_cur_time << std::endl;
        // clang-format on
    } else {
        m_cur_time = m_new_time;
    }

    return continue_sim;
}

void SimTime::set_current_cfl(amrex::Real cfl_unit_time)
{
    amrex::Real dt_new = 2.0 * m_max_cfl / cfl_unit_time;

    // Restrict timestep during initialization phase
    if (m_is_init) dt_new *= m_init_shrink;

    // Limit timestep growth to 10% per timestep
    if (m_dt[0] > 0.0) dt_new = amrex::min(dt_new, 1.1 * m_dt[0]);

    // Don't overshoot stop time
    if ((m_stop_time > 0.0) && ((m_cur_time + dt_new) > m_stop_time))
        dt_new = m_stop_time - m_cur_time;

    if (m_adaptive) {
        m_dt[0] = dt_new;

        if(m_is_init && m_initial_dt > 0.0) m_dt[0] = m_initial_dt;

        if (!m_is_init) m_new_time = m_cur_time + m_dt[0];

    } else {
        // If user has specified fixed DT then issue a warning if the timestep
        // is larger than the deltaT determined from max. CFL considerations.
        if ((dt_new < m_fixed_dt) && !m_is_init) {
            amrex::Print()
                << "WARNING: fixed_dt does not satisfy CFL condition.\n"
                << "deltaT (max. CFL) = " << dt_new
                << " fixed_dt = " << m_fixed_dt << std::endl;
        }
        // Ensure that we use user-specified dt. Checkpoint restart might have
        // overridden this
        m_dt[0] = m_fixed_dt;
    }

    m_current_cfl = 0.5 * cfl_unit_time * m_dt[0];
    if (m_verbose >= 0)
        amrex::Print() << "CFL: " << m_current_cfl << " dt: " << m_dt[0]
                       << std::endl;
}

bool SimTime::continue_simulation()
{
    constexpr double eps = 1.0e-12;
    bool stop_simulation = false;

    if (m_stop_time_index == 0) return stop_simulation;

    if ((m_stop_time > 0.0) && ((m_new_time + eps) >= m_stop_time))
        return stop_simulation;

    if ((m_stop_time_index > 0) && (m_time_index >= m_stop_time_index))
        return stop_simulation;

    return !(stop_simulation);
}

bool SimTime::do_regrid()
{
    return (
        (m_regrid_interval > 0) && (m_time_index > 0) &&
        (m_time_index % m_regrid_interval == 0));
}

bool SimTime::write_plot_file()
{
    return ((m_plt_interval > 0) && (m_time_index % m_plt_interval == 0));
}

bool SimTime::write_checkpoint()
{
    return ((m_chkpt_interval > 0) && (m_time_index % m_chkpt_interval == 0));
}

bool SimTime::write_last_plot_file()
{
    return ((m_plt_interval > 0) && (m_time_index % m_plt_interval != 0));
}

bool SimTime::write_last_checkpoint()
{
    return ((m_chkpt_interval > 0) && (m_time_index % m_chkpt_interval != 0));
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
