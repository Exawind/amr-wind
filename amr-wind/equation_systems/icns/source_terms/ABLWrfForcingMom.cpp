#include "amr-wind/equation_systems/icns/source_terms/ABLForcing.H"
#include "amr-wind/utilities/PlaneAveraging.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABL.H"


#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"

namespace amr_wind{
namespace pde{
namespace icns{

ABLWrfForcingMom::ABLWrfForcingMom(const CFDSim& sim) : m_time(sim.time())
{

  const auto& abl = sim.physics_manager().get<amr_wind::ABL>();
  amrex::ParmParse pp_abl(identifier());

  pp_abl.get("WRF force file", m_wrf_file);

}

ABLWrfForcingMom::~ABLWrfForcingMom = default;

void ABLWrfForcingMom::read_forcing_file()
{
  auto ncf = ncutils::NCFile::open_par(
      m_filename, NC_NOWRITE | NC_NETCDF4 | NC_MPIIO,
      amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

  m_nheight  = ncf.dim("nheight").len();
  m_ntime = ncf.dim("ntime").len();

  m_wrf_height.reize(m_nheight);
  m_wrf_time.resize(m_ntime);

  ncf.var("heights").get(m_wrf_height.data());
  ncf.var("times").get(m_wrf_time.data());

  m_wrf_mom

}
                                                        
} // namespace icns
} // namespace pde
} // namespace amr_wind
