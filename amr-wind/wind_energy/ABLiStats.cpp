#include "amr-wind/wind_energy/ABLiStats.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {

ABLiStats::ABLiStats(CFDSim& sim, const ABLWallFunction& abl_wall_func)
    : m_sim(sim), m_abl_wall_func(abl_wall_func)
{}

ABLiStats::~ABLiStats() = default;

void ABLiStats::initialize()
{
    BL_PROFILE("amr-wind::ABLiStats::initialize");

    // Labels for the different sampler types
    amrex::Vector<std::string> labels;
    // Fields to be sampled - requested by user
    amrex::Vector<std::string> field_names;

    {
        amrex::ParmParse pp("ABL");
        pp.query("stats_output_frequency", m_out_freq);
        pp.query("stats_output_format", m_out_fmt);
    }

    if (m_out_fmt == "netcdf") prepare_netcdf_file();
}

void ABLiStats::post_advance_work()
{
    BL_PROFILE("amr-wind::ABLiStats::post_advance_work");
    const auto& time = m_sim.time();
    const int tidx = time.time_index();
    // Skip processing if it is not an output timestep
    if (!(tidx % m_out_freq == 0)) return;

    compute_zi();
    process_output();
}

void ABLiStats::compute_zi() {

    //TODO: Implement this function
    m_zi = 0.0;
}

void ABLiStats::process_output()
{
    if (m_out_fmt == "ascii") {
        write_ascii();
    } else if (m_out_fmt == "netcdf") {
        write_netcdf();
    } else {
        amrex::Abort("ABLiStats: Invalid output format encountered");
    }
}

void ABLiStats::write_ascii()
{
    BL_PROFILE("amr-wind::ABLiStats::write_ascii");
    amrex::Print() << "WARNING: ABLiStats: ASCII output will impact performance"
                   << std::endl;

    const std::string istat_dir = "i_stats";
    const std::string sname = amrex::Concatenate("i_stats", m_sim.time().time_index());

    if (!amrex::UtilCreateDirectory(istat_dir, 0755)) {
        amrex::CreateDirectoryFailed(istat_dir);
    }
    const std::string fname = istat_dir + "/" + sname + ".txt";
    //TODO: Do something to write the data to file here
}

void ABLiStats::prepare_netcdf_file()
{
#ifdef AMR_WIND_USE_NETCDF

    const std::string istat_dir = "i_stats";
    const std::string sname = amrex::Concatenate("i_stats", m_sim.time().time_index());
    if (!amrex::UtilCreateDirectory(istat_dir, 0755)) {
        amrex::CreateDirectoryFailed(istat_dir);
    }
    m_ncfile_name = istat_dir + "/" + sname + ".nc";

    // Only I/O processor handles NetCDF generation
    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    auto ncf = ncutils::NCFile::create(m_ncfile_name, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "num_time_steps";
    const std::string npart_name = "num_points";
    const std::vector<std::string> two_dim{nt_name, npart_name};
    ncf.enter_def_mode();
    ncf.put_attr("title", "AMR-Wind data sampling output");
    ncf.put_attr("version", ioutils::amr_wind_version());
    ncf.put_attr("created_on", ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim("ndim", AMREX_SPACEDIM);
    ncf.def_var("time", NC_DOUBLE, {nt_name});
    ncf.def_var("ustar", NC_DOUBLE, {nt_name});
    ncf.def_var("wstar", NC_DOUBLE, {nt_name});
    ncf.def_var("L", NC_DOUBLE, {nt_name});
    ncf.def_var("zi", NC_DOUBLE, {nt_name});
    ncf.exit_def_mode();

#else
    amrex::Abort(
        "NetCDF support was not enabled during build time. Please recompile or "
        "use native format");
#endif
}

void ABLiStats::write_netcdf()
{
#ifdef AMR_WIND_USE_NETCDF

    if (!amrex::ParallelDescriptor::IOProcessor()) return;
    auto ncf = ncutils::NCFile::open(m_ncfile_name, NC_WRITE);
    const std::string nt_name = "num_time_steps";
    // Index of the next timestep
    const size_t nt = ncf.dim(nt_name).len();
    {
        auto time = m_sim.time().new_time();
        ncf.var("time").put(&time, {nt}, {1});
        auto ustar = m_abl_wall_func.utau();
        ncf.var("ustar").put(&ustar, {nt}, {1});
        double wstar = 0.0;
        ncf.var("wstar").put(&wstar, {nt}, {1});
        double L = 1e18;
        ncf.var("L").put(&L, {nt}, {1});
        ncf.var("zi").put(&m_zi, {nt}, {1});
    }
    ncf.close();
#endif
}

}
