#include "amr-wind/wind_energy/actuator/turbine/external/ExtTurbIface.H"
#include "amr-wind/wind_energy/actuator/turbine/kynema/kynema_types.H"
#include "amr-wind/wind_energy/actuator/turbine/kynema/kynema_wrapper.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_FileSystem.H"

#include <algorithm>
#include <cmath>

namespace ext_turb {

template <>
ExtTurbIface<KynemaTurbine, KynemaSolverData>::~ExtTurbIface()
{
    //! Do deallocation if necessary
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::parse_inputs(
    const amr_wind::CFDSim& sim, const std::string& inp_name)
{
    amrex::ParmParse pp(inp_name);

    const auto& time = sim.time();
    //! Check that the user has not enabled adaptive timestepping
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !time.adaptive_timestep(),
        "Adaptive time-stepping not supported when Kynema is enabled");

    m_dt_cfd = time.delta_t();

    //! Put input-reading steps here
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::allocate_ext_turbines()
{
    BL_PROFILE("amr-wind::KynemaIface::allocate_turbines");
    int nturbines = static_cast<int>(m_turbine_data.size());
    //!! Allocate stuff in Kynema
    m_is_initialized = true;
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::init_solution(
    const int local_id)
{
    BL_PROFILE("amr-wind::KynemaIface::init_solution");
    AMREX_ALWAYS_ASSERT(local_id < static_cast<int>(m_turbine_data.size()));
    AMREX_ALWAYS_ASSERT(m_is_initialized);

    auto& fi = *m_turbine_data[local_id];
    //!! Set up initial solution
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::get_hub_stats(
    const int local_id)
{
    BL_PROFILE("amr-wind::KynemaIface::get_hub_stats");

    auto& fi = *m_turbine_data[local_id];
}

#ifdef AMR_WIND_USE_KYNEMA

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::prepare_netcdf_file(
    KynemaTurbine& fi)
{
#ifdef AMR_WIND_USE_NETCDF
    BL_PROFILE("amr-wind::KynemaIface::prepare_netcdf_file");
    if (!amrex::UtilCreateDirectory(m_solver_data.m_output_dir, 0755)) {
        amrex::CreateDirectoryFailed(m_solver_data.m_output_dir);
    }

    const std::string fname =
        m_solver_data.m_output_dir + "/" + fi.tlabel + ".nc";

    // Don't overwrite existing
    if (amrex::FileSystem::Exists(fname)) {
        return;
    }

    auto ncf = ncutils::NCFile::create(fname, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "num_time_steps";
    const std::string np_name = "num_vel_points";
    ncf.enter_def_mode();
    ncf.put_attr("title", "AMR-Wind Kynema velocity data");
    ncf.put_attr("AMR-Wind_version", amr_wind::ioutils::amr_wind_version());
    ncf.put_attr("created_on", amr_wind::ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim(np_name, fi.from_cfd.u_Len);
    ncf.def_dim("ndim", AMREX_SPACEDIM);
    ncf.def_var("time", NC_FLOAT, {nt_name});
    ncf.def_var("xco", NC_FLOAT, {np_name});
    ncf.def_var("yco", NC_FLOAT, {np_name});
    ncf.def_var("zco", NC_FLOAT, {np_name});
    ncf.def_var("uvel", NC_FLOAT, {nt_name, np_name});
    ncf.def_var("vvel", NC_FLOAT, {nt_name, np_name});
    ncf.def_var("wvel", NC_FLOAT, {nt_name, np_name});
    ncf.exit_def_mode();

    {
        const auto npts = static_cast<size_t>(fi.from_cfd.u_Len);
        auto xco = ncf.var("xco");
        xco.put(fi.position_at_vel(0), {0}, {npts});
        auto yco = ncf.var("yco");
        yco.put(fi.position_at_vel(1), {0}, {npts});
        auto zco = ncf.var("zco");
        zco.put(fi.position_at_vel(2), {0}, {npts});
    }
#else
    amrex::ignore_unused(fi);
    amrex::OutStream()
        << "WARNING: KynemaIface: NetCDF support was not enabled during compile "
           "time. KynemaIface cannot support restart."
        << std::endl;
#endif
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::write_velocity_data(
    const KynemaTurbine& fi)
{
#ifdef AMR_WIND_USE_NETCDF
    BL_PROFILE("amr-wind::KynemaIface::write_velocity_data");
    const std::string fname =
        m_solver_data.m_output_dir + "/" + fi.tlabel + ".nc";
    auto ncf = ncutils::NCFile::open(fname, NC_WRITE);
    const std::string nt_name = "num_time_steps";
    const size_t nt = ncf.dim(nt_name).len();
    const auto npts = static_cast<size_t>(fi.from_cfd.u_Len);

    const double time = fi.time_index * fi.dt_ext;
    ncf.var("time").put(&time, {nt}, {1});
    const auto& uu = fi.from_cfd;
    ncf.var("uvel").put(uu.u, {nt, 0}, {1, npts});
    ncf.var("vvel").put(uu.v, {nt, 0}, {1, npts});
    ncf.var("wvel").put(uu.w, {nt, 0}, {1, npts});
#else
    amrex::ignore_unused(fi);
#endif
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::read_velocity_data(
    KynemaTurbine& fi, const ncutils::NCFile& ncf, const size_t tid)
{
#ifdef AMR_WIND_USE_NETCDF
    const auto nt = static_cast<size_t>(tid);
    const auto npts = static_cast<size_t>(fi.from_cfd.u_Len);

    auto& uu = fi.from_cfd;
    ncf.var("uvel").get(uu.u, {nt, 0}, {1, npts});
    ncf.var("vvel").get(uu.v, {nt, 0}, {1, npts});
    ncf.var("wvel").get(uu.w, {nt, 0}, {1, npts});
#else
    amrex::ignore_unused(fi);
    amrex::Abort(
        "KynemaIface::read_velocity_data: AMR-Wind was not compiled with NetCDF "
        "support");
#endif
}

#else

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::prepare_netcdf_file(
    KynemaTurbine& /*unused*/)
{}
template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::write_velocity_data(
    const KynemaTurbine& /*unused*/)
{}
template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::read_velocity_data(
    KynemaTurbine& /*unused*/,
    const ncutils::NCFile& /*unused*/,
    const size_t /*unused*/)
{}

#endif

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::do_turbine_step(int& tid)
{
    // individual turbine step
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::write_turbine_checkpoint(
    int& tid)
{
    // write checkpoint
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::ext_init_turbine(
    KynemaTurbine& fi)
{
    BL_PROFILE("amr-wind::KynemaIface::init_turbine");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        amrex::FileSystem::Exists(fi.input_file),
        "KynemaIface: Cannot find Kynema input file: " + fi.input_file);

    // start up turbine

    // Determine the number of substeps for Kynema per CFD timestep
    fi.num_substeps = static_cast<int>(std::floor(fi.dt_cfd / fi.dt_ext));

    AMREX_ALWAYS_ASSERT(fi.num_substeps > 0);
    // Check that the time step sizes are consistent and Kynema advances at an
    // integral multiple of CFD timestep
    double dt_err =
        fi.dt_cfd / (static_cast<double>(fi.num_substeps) * fi.dt_ext) - 1.0;
    if (dt_err > 1.0e-12) {
        amrex::Abort(
            "KynemaIFace: Kynema timestep is not an integral "
            "multiple of CFD timestep");
    }
}

// cppcheck-suppress constParameterReference
// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::ext_replay_turbine(
    KynemaTurbine& fi)
{

    // Do we even do this???
    
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::ext_restart_turbine(
    KynemaTurbine& fi)
{
    BL_PROFILE("amr-wind::KynemaIface::restart_turbine");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        amrex::FileSystem::Exists(fi.checkpoint_file + ".chkp"),
        "KynemaIface: Cannot find Kynema checkpoint file: " +
            fi.checkpoint_file);

    // Determine the number of substeps for Kynema per CFD timestep
    fi.num_substeps = static_cast<int>(std::floor(fi.dt_cfd / fi.dt_ext));

    AMREX_ALWAYS_ASSERT(fi.num_substeps > 0);
    // Check that the time step sizes are consistent and Kynema advances at an
    // integral multiple of CFD timestep
    double dt_err =
        fi.dt_cfd / (static_cast<double>(fi.num_substeps) * fi.dt_ext) - 1.0;
    if (dt_err > 1.0e-4) {
        amrex::Abort(
            "KynemaIFace: Kynema timestep is not an integral "
            "multiple of CFD timestep");
    }
}

template class ExtTurbIface<KynemaTurbine, KynemaSolverData>;

} // namespace ext_turb
