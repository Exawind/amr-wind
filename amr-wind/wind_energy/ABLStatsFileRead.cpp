#include "amr-wind/wind_energy/ABLStatsFileRead.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_Print.H"
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

namespace amr_wind {

ABLReadStats::ABLReadStats(const std::string filestats)
    : m_stat_filename(filestats)
{

    auto ncf = ncutils::NCFile::open_par(
        m_stat_filename, NC_NOWRITE | NC_NETCDF4 | NC_MPIIO,
        amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

    m_stats_nt_steps = ncf.dim("time").len();

    m_stats_time.resize(m_stats_nt_steps);

    m_stats_ustar.resize(m_stats_nt_steps);

    ncf.var("time").get(m_stats_time.data());

    ncf.var("ustar").get(m_stats_ustar.data());

    auto grp = ncf.group("mean_profiles");

    m_stats_nlevels = grp.dim("nlevels").len();

    m_stats_theta.resize(m_stats_nt_steps * m_stats_nlevels);

    grp.var("theta").get(m_stats_theta.data());

    ncf.close();
}

const amrex::Vector<amrex::Real>& ABLReadStats::stats_theta() const
{
    return m_stats_theta;
}

const amrex::Vector<amrex::Real>& ABLReadStats::stats_time() const
{
    return m_stats_time;
}

const amrex::Vector<amrex::Real>& ABLReadStats::stats_ustar() const
{
    return m_stats_ustar;
}

} // namespace amr_wind
