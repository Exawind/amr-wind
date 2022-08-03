#include "amr-wind/wind_energy/ABLMesoscaleInput.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_Print.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

ABLWRFfile::ABLWRFfile(const std::string wrfFile) : m_wrf_filename(wrfFile)
{
#ifdef AMR_WIND_USE_NETCDF
    auto ncf = ncutils::NCFile::open_par(
        m_wrf_filename, NC_NOWRITE | NC_NETCDF4 | NC_MPIIO,
        amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

    m_wrf_nheight = ncf.has_dim("nheight") ? ncf.dim("nheight").len() : 0;
    m_wrf_ntime = ncf.dim("ntime").len();
    amrex::Print() << "Loading " << m_wrf_filename << " : " << m_wrf_ntime
                   << " times, " << m_wrf_nheight << " heights" << std::endl;

    m_wrf_height.resize(m_wrf_nheight);
    m_wrf_time.resize(m_wrf_ntime);

    if (m_wrf_nheight > 0) ncf.var("heights").get(m_wrf_height.data());
    ncf.var("times").get(m_wrf_time.data());

    m_wrf_u.resize(m_wrf_nheight * m_wrf_ntime);
    m_wrf_v.resize(m_wrf_nheight * m_wrf_ntime);
    m_wrf_temp.resize(m_wrf_nheight * m_wrf_ntime);
    m_wrf_tflux.resize(m_wrf_ntime);

    if (m_wrf_nheight > 0) {
        ncf.var("wrf_momentum_u").get(m_wrf_u.data());
        ncf.var("wrf_momentum_v").get(m_wrf_v.data());
        ncf.var("wrf_temperature").get(m_wrf_temp.data());
    } else {
        amrex::Print() << "No height dimension in netcdf input file; no "
                          "forcing profiles read."
                       << std::endl;
    }
    ncf.var("wrf_tflux").get(m_wrf_tflux.data());

    // ***FIXME***
    // MUST COMMENT THIS LINE OUT (resize cmd) to consistently fix problem:
    //
    // m_wrf_transition_height.resize(m_wrf_ntime);
    //
    // if (ncf.has_var("transition_height")) {
    //    amrex::Print() << "found transition_height in WRFforcing file" <<
    //    std::endl;
    //    ncf.var("transition_height").get(m_wrf_transition_height.data());
    //}

#else
    amrex::Abort(
        "NetCDF support was not enabled during build time. Please recompile");
#endif

    amrex::ParmParse pp("ABL");
    pp.query("WRF_tendency_forcing", m_abl_wrf_tendency);
}

const amrex::Vector<amrex::Real>& ABLWRFfile::wrf_heights() const
{
    return m_wrf_height;
}

const amrex::Vector<amrex::Real>& ABLWRFfile::wrf_times() const
{
    return m_wrf_time;
}

const amrex::Vector<amrex::Real>& ABLWRFfile::wrf_u() const { return m_wrf_u; }

const amrex::Vector<amrex::Real>& ABLWRFfile::wrf_v() const { return m_wrf_v; }

const amrex::Vector<amrex::Real>& ABLWRFfile::wrf_temp() const
{
    return m_wrf_temp;
}

const amrex::Vector<amrex::Real>& ABLWRFfile::wrf_tflux() const
{
    return m_wrf_tflux;
}

// ***FIXME***
// const amrex::Vector<amrex::Real>& ABLWRFfile::wrf_transition_height() const {
// return m_wrf_transition_height; }

bool ABLWRFfile::is_wrf_tendency_forcing() const { return m_abl_wrf_tendency; }

int ABLWRFfile::nheights() const { return m_wrf_nheight; }
int ABLWRFfile::times() const { return m_wrf_ntime; }

} // namespace amr_wind
