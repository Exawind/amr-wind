#include "amr-wind/wind_energy/ABLMesoscaleInput.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_Print.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

// cppcheck-suppress uninitMemberVar
ABLMesoscaleInput::ABLMesoscaleInput(
    const std::string& ncfile, const std::string& var_prefix)
    : m_filename(ncfile)
{
#ifdef AMR_WIND_USE_NETCDF
    auto ncf = ncutils::NCFile::open_par(
        m_filename, NC_NOWRITE | NC_NETCDF4 | NC_MPIIO,
        amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

    m_nheight = ncf.has_dim("nheight") ? ncf.dim("nheight").len() : 0;
    m_ntime = ncf.dim("ntime").len();
    amrex::Print() << "Loading " << m_filename << " : " << m_ntime << " times, "
                   << m_nheight << " heights" << std::endl;

    m_height.resize(m_nheight);
    m_time.resize(m_ntime);

    if (m_nheight > 0) ncf.var("heights").get(m_height.data());
    ncf.var("times").get(m_time.data());

    m_u.resize(m_nheight * m_ntime);
    m_v.resize(m_nheight * m_ntime);
    m_temp.resize(m_nheight * m_ntime);
    m_tflux.resize(m_ntime);

    if (m_nheight > 0) {
        ncf.var(var_prefix + "momentum_u").get(m_u.data());
        ncf.var(var_prefix + "momentum_v").get(m_v.data());
        ncf.var(var_prefix + "temperature").get(m_temp.data());
    } else {
        amrex::Print() << "No height dimension in netcdf input file; no "
                          "forcing profiles read."
                       << std::endl;
    }
    ncf.var(var_prefix + "tflux").get(m_tflux.data());

    // ***FIXME***
    // MUST COMMENT THIS LINE OUT (resize cmd) to consistently fix problem:
    //
    // m_transition_height.resize(m_ntime);
    //
    // if (ncf.has_var("transition_height")) {
    //    amrex::Print() << "found transition_height in ABLMesoscaleInput file"
    //                   << std::endl;
    //    ncf.var("transition_height").get(m_transition_height.data());
    //}

#else
    amrex::Abort(
        "NetCDF support was not enabled during build time. Please recompile");
#endif

    amrex::ParmParse pp("ABL");
    pp.query("tendency_forcing", m_abl_tendency);
}

} // namespace amr_wind
