#ifndef FAST_TYPES_H
#define FAST_TYPES_H

#include "amr-wind/wind_energy/actuator/turbine/fast/fast_wrapper.H"

#include <string>

namespace exw_fast {

enum class SimMode : int {
    init = 0, ///< Clean start
    replay,   ///< Replay using velocities stored in file
    restart   ///< Restart using FAST checkpoint files
};

/** Representation of a turbine for exchanging data with FAST
 */
struct FastTurbine
{
    static constexpr int ndim = 3;

    //! Unique string identifier for this turbine
    std::string tlabel;

    //! Local ID for this turbine (provided by FAST interface)
    int tid_local;

    //! Global ID for this turbine (set by Turbine instance)
    int tid_global;

    //! Number of actuator points per blade
    int num_pts_blade;

    //! Number of actuator points for tower
    int num_pts_tower;

    //! Position of tower base in global coordinate system
    float base_pos[ndim];

    SimMode sim_mode{SimMode::init};

    //! OpenFAST input file
    std::string input_file;

    //! Checkpoint file name
    std::string checkpoint_file;

    //! Number of blades
    int num_blades;

    //! Total number of elements along the blade
    int num_blade_elem;

    //! Start time for this turbine
    double start_time{0.0};

    //! End time for this turbine
    double stop_time;

    //! Timestep for CFD
    double dt_cfd;

    //! Timestep for FAST using per-turbine to check input files
    double dt_fast;

    //! Number of sub-steps of fast per CFD timestep
    int num_substeps;

    //! Time step index for FAST
    int time_index{0};

    //! Does FAST need solution0
    bool is_solution0{true};

    //! Checkpoint interval for FAST
    int chkpt_interval;

    // Data structures that are used to exchange between fast/cfd

    exw_fast::OpFM_InputType to_cfd;
    exw_fast::OpFM_OutputType from_cfd;

#ifdef AMR_WIND_FAST_USE_SCDX
    exw_fast::SC_DX_InputType to_sc;
    exw_fast::SC_DX_OutputType from_sc;
#else
    exw_fast::SC_InputType to_sc;
    exw_fast::SC_OutputType from_sc;
#endif
};

} // namespace exw_fast

#endif /* FAST_TYPES_H */
