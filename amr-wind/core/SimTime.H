#ifndef SIMTIME_H
#define SIMTIME_H

#include "AMReX_REAL.H"
#include "AMReX_GpuQualifiers.H"
#include "AMReX_Extension.H"

namespace amr_wind {

/** Time manager for simulations
 *  \ingroup core
 */
class SimTime
{
public:
    static constexpr int max_time_states = 3;

    SimTime() = default;

    ~SimTime() = default;

    /** Advance to a new time state and update internal time states
     *
     *  Return a boolean flag indicating if the simulation should continue.
     */
    bool new_timestep();

    /** Return true if simulation should continue
     */
    bool continue_simulation();

    /** Return true if mesh refinement should be performed at this timestep
     */
    bool do_regrid();

    /** Return true if plot files should be written at this timestep
     */
    bool write_plot_file();

    /** Return true if checkpoint files should be written at this timestep
     */
    bool write_checkpoint();

    /** Return true if plot files should be written at this timestep
     */
    bool write_last_plot_file();

    /** Return true if checkpoint files should be written at this timestep
     */
    bool write_last_checkpoint();

    /** Set current CFL and update timestep based on CFL
     *
     *  @param cfl_unit_time `(U / dx)`, i.e., CFL assuming unit dt
     */
    void set_current_cfl(amrex::Real cfl_unit_time);

    /** Set start time and index based on checkpoint/restart file
     *
     *  @param tidx Time index from checkpoint file
     *  @param time Time read from checkpoint file
     */
    void set_restart_time(int tidx, amrex::Real time);

    AMREX_FORCE_INLINE
    amrex::Real deltaT() const { return m_dt[0]; }

    AMREX_FORCE_INLINE
    amrex::Real& deltaT() { return m_dt[0]; }

    AMREX_FORCE_INLINE
    amrex::Real deltaTNm1() const { return m_dt[1]; }

    AMREX_FORCE_INLINE
    amrex::Real deltaTNm2() const { return m_dt[2]; }

    AMREX_FORCE_INLINE
    amrex::Real& deltaTNm1() { return m_dt[1]; }

    AMREX_FORCE_INLINE
    amrex::Real& deltaTNm2() { return m_dt[2]; }

    AMREX_FORCE_INLINE
    amrex::Real current_time() const { return m_cur_time; }

    AMREX_FORCE_INLINE
    amrex::Real& current_time() { return m_cur_time; }

    AMREX_FORCE_INLINE
    amrex::Real new_time() const { return m_new_time; }

    AMREX_FORCE_INLINE
    amrex::Real max_cfl() const { return m_max_cfl; }

    AMREX_FORCE_INLINE
    int time_index() const { return m_time_index; }

    AMREX_FORCE_INLINE
    int& time_index() { return m_time_index; }

    AMREX_FORCE_INLINE
    bool adaptive_timestep() const { return m_adaptive; }

    AMREX_FORCE_INLINE
    bool use_force_cfl() const { return m_use_force_cfl; }

    //! Read user defined options from input file
    void parse_parameters();

private:
    //! Timestep sizes
    amrex::Real m_dt[max_time_states];

    //! Current simulation time
    amrex::Real m_cur_time{0.0};

    //! Time where a new timestep will advance to
    amrex::Real m_new_time{0.0};

    //! Starting simulation time (to track restart)
    amrex::Real m_start_time{0.0};

    //! Max time for simulation
    amrex::Real m_stop_time{-1.0};

    //! Maximum CFL constraint
    amrex::Real m_max_cfl{0.5};

    //! Fixed deltaT indicated by user
    amrex::Real m_fixed_dt{-1.0};

    //! Initial deltaT indicated by user
    amrex::Real m_initial_dt{-1.0};

    //! Current CFL
    amrex::Real m_current_cfl{0.0};

    //! Initial reduction in timestep size for startup
    amrex::Real m_init_shrink{0.1};

    //! Counter for the number of timesteps since start of simulation
    int m_time_index{0};

    //! Initial starting time index
    int m_start_time_index{0};

    //! Maximum timesteps for simulation
    int m_stop_time_index{-1};

    //! Time interval for plot file output
    int m_plt_interval{-1};

    //! Time interval for writing checkpoint/restart files
    int m_chkpt_interval{-1};

    //! Time interval for regridding
    int m_regrid_interval{-1};

    //! Verbosity
    int m_verbose{0};

    //! Flag indicating whether timestep is adaptive
    bool m_adaptive{false};

    //! Flag indicating if this is initialization
    bool m_is_init{true};

    //! Flag indicating if forcing should be included in CFL calculation
    bool m_use_force_cfl{true};

};

} // namespace amr_wind

#endif /* SIMTIME_H */
