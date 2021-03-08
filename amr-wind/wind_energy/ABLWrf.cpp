#include "amr-wind/wind_energy/ABLWrf.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "AMReX_ParallelDescriptor.H"

namespace amr_wind {

ABLWRFfile::ABLWRFfile(const std::string filewrf)
    : m_wrf_filename(filewrf)
{

    auto ncf = ncutils::NCFile::open_par(
        m_wrf_filename, NC_NOWRITE | NC_NETCDF4 | NC_MPIIO,
        amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

    m_nheight = ncf.dim("nheight").len();
    m_ntime = ncf.dim("ntime").len();

    m_wrf_height.resize(m_nheight);
    m_wrf_time.resize(m_ntime);

    ncf.var("heights").get(m_wrf_height.data());
    ncf.var("times").get(m_wrf_time.data());

    m_wrf_u.resize(m_nheight * m_ntime);
    m_wrf_v.resize(m_nheight * m_ntime);
    m_wrf_temp.resize(m_nheight * m_ntime);
    m_wrf_tflux.resize(m_ntime);

    ncf.var("wrf_momentum_u").get(m_wrf_u.data());
    ncf.var("wrf_momentum_v").get(m_wrf_v.data());
    ncf.var("wrf_temperature").get(m_wrf_temp.data());
    ncf.var("wrf_tflux").get(m_wrf_tflux.data());
}

    const amrex::Vector<amrex::Real>& wrf_heights() const
    {
        return m_wrf_height;
    }

    const amrex::Vector<amrex::Real>& wrf_times() const { return m_wrf_time; }

    const amrex::Vector<amrex::Real>& wrf_u() const { return m_wrf_u; }

    const amrex::Vector<amrex::Real>& wrf_v() const { return m_wrf_v; }

    const amrex::Vector<amrex::Real>& wrf_temp() const { return m_wrf_temp; }

    const amrex::Vector<amrex::Real>& wrf_tflux() const { return m_wrf_tflux; }

    int nheights() const { return m_nheight; }
    int times() const { return m_ntime; }


} // namespace amr_wind
