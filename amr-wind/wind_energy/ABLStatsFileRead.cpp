#include "amr-wind/wind_energy/ABLStatsFileRead.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_Print.H"
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <cstddef>

namespace amr_wind {

namespace {

//! Return closest index (from lower) of value in vector
AMREX_FORCE_INLINE int
closest_index(const amrex::Vector<amrex::Real>& vec, const amrex::Real value)
{
    auto const it = std::upper_bound(vec.begin(), vec.end(), value);
    AMREX_ALWAYS_ASSERT(it != vec.end());

    const int idx = std::distance(vec.begin(), it);
    return std::max(idx - 1, 0);
}
} // namespace

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
    m_stats_theta_1D.resize(m_stats_nlevels);

    m_stats_u.resize(m_stats_nt_steps*m_stats_nlevels);
    m_stats_u_1D.resize(m_stats_nlevels);
    
    m_stats_v.resize(m_stats_nt_steps*m_stats_nlevels);
    m_stats_v_1D.resize(m_stats_nlevels);
    
    m_stats_vmag.resize(m_stats_nt_steps*m_stats_nlevels);
    m_stats_vmag_1D.resize(m_stats_nlevels);

    grp.var("theta").get(m_stats_theta.data());
    grp.var("hvelmag").get(m_stats_vmag.data());
    grp.var("u").get(m_stats_u.data());
    grp.var("v").get(m_stats_v.data()); 

    ncf.close();
}

const amrex::Vector<amrex::Real>& ABLReadStats::stats_theta() const
{
    return m_stats_theta_1D;
}

const amrex::Vector<amrex::Real>& ABLReadStats::stats_time() const
{
    return m_stats_time;
}

const amrex::Vector<amrex::Real>& ABLReadStats::stats_ustar() const
{
    return m_stats_ustar;
}

amrex::Real ABLReadStats::interpUstarTime(const amrex::Real timeSim)
{

    // First the index in time
    int idx_time;
    idx_time = closest_index(m_stats_time, timeSim);

    amrex::Array<amrex::Real, 2> coeff_interp{{0.0, 0.0}};

    amrex::Real denom = m_stats_time[idx_time + 1] - m_stats_time[idx_time];

    coeff_interp[0] = (m_stats_time[idx_time + 1] - timeSim) / denom;
    coeff_interp[1] = 1.0 - coeff_interp[0];

    amrex::Real interpUstar;

    interpUstar = coeff_interp[0] * m_stats_time[idx_time] +
                  coeff_interp[1] * m_stats_time[idx_time + 1];

    return interpUstar;
}
void ABLReadStats::interpThetaTime(const amrex::Real timeSim, amrex::Array<amrex::Real, 4>& wall_func_aux)
{

    // First the index in time
    int idx_time;
    idx_time = closest_index(m_stats_time, timeSim);

    amrex::Array<amrex::Real, 2> coeff_interp{{0.0, 0.0}};

    amrex::Real denom = m_stats_time[idx_time + 1] - m_stats_time[idx_time];

    coeff_interp[0] = (m_stats_time[idx_time + 1] - timeSim) / denom;
    coeff_interp[1] = 1.0 - coeff_interp[0];

    for (size_t i = 0; i < m_stats_nlevels; i++) {
        size_t lt = idx_time * m_stats_nlevels + i;
        size_t rt = (idx_time + 1) * m_stats_nlevels + i;

        m_stats_theta_1D[i] = coeff_interp[0] * m_stats_theta[lt] +
                              coeff_interp[1] * m_stats_theta[rt];

        m_stats_u_1D[i] =
            coeff_interp[0] * m_stats_u[lt] + coeff_interp[1] * m_stats_u[rt];

        m_stats_v_1D[i] =
            coeff_interp[0] * m_stats_v[lt] + coeff_interp[1] * m_stats_v[rt];

        m_stats_vmag_1D[i] = coeff_interp[0] * m_stats_vmag[lt] +
                             coeff_interp[1] * m_stats_vmag[rt];
    }
    wall_func_aux[0] = m_stats_u_1D[0];
    wall_func_aux[1] = m_stats_v_1D[0];
    wall_func_aux[2] = m_stats_vmag_1D[0];
    wall_func_aux[3] = m_stats_theta_1D[0];
    
}

} // namespace amr_wind
