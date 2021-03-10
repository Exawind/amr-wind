#include "amr-wind/equation_systems/icns/source_terms/ABLWrfForcingTemp.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Print.H"

namespace amr_wind{
namespace pde{
namespace icns{

ABLWrfForcingTemp::ABLWrfForcingTemp(const CFDSim& sim)
    : m_time(sim.time()), m_mesh(sim.mesh())
{
  const auto& abl = sim.physics_manager().get<amr_wind::ABL>();
  abl.register_mean_wrf_temp_forcing(this);

  mean_temperature_init(abl.abl_statistics().theta_profile, abl.abl_wrf_file);
}

ABLWrfForcingTemp::~ABLWrfForcingTemp() = default;

void ABLWrfForcingTemp::mean_temperature_init(const FieldPlaneAveraging& tavg, const ABLWRFfile& wrfFile)
{
  m_axis = vavg.axis();
  // The implementation depends the assumption that the ABL statistics class
  // computes statistics at the cell-centeres only on level 0. If this
  // assumption changes in future, the implementation will break... so put in
  // a check here to catch this.
  AMREX_ALWAYS_ASSERT(
      m_mesh.Geom(0).Domain().length(m_axis) ==
      static_cast<int>(vavg.line_centroids().size()));

  m_theta_ht.resize(tavg.line_centroids().size());
  m_theta_vals.resize(tavg.ncell_line());

  amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, tavg.line_centroids().begin(),
      tavg.line_centroids().end(), m_theta_ht.begin());
  
  m_wrf_theta_vals.resize(wrfFile.nheights());
  m_wrf_ht.resize(wrfFile.nheights());

  amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, wrfFile.wrf_heights().begin(),
      wrfFile.wrf_heights().end(), m_wrf_ht.begin());
    
}

void ABLWrfForcingTemp::mean_temperature_heights(const FieldPlaneAveraging& tavg, std::unique_ptr<ABLWRFfile>& wrfFile)
{

  amrex::Real currtime;
  currtime = m_time.current_time();

  // First the index in time
  m_idx_time = closest_index(wrfFile->wrf_times(), currtime);

  amrex::Array<amrex::Real, 2> coeff_interp{{0.0, 0.0}};

   amrex::Real denom  = wrfFile->wrf_times()[m_idx_time+1] -
       wrfFile->wrf_times()[m_idx_time];
  
  coeff_interp[0] =  (wrfFile->wrf_times()[m_idx_time+1] - currtime) / denom; 
  coeff_interp[1] = 1.0 - coeff_interp[0];

  int num_wrf_ht = wrfFile->nheights();
  
  amrex::Vector<amrex::Real> wrfInterptheta(num_wrf_ht);

  for (int i = 0 ; i < num_wrf_ht; i++) {
    int lt = m_idx_time*num_wrf_ht + i;
    int rt = (m_idx_time+1)*num_wrf_ht + i;
    
     wrfInterptheta[i] = coeff_interp[0]*wrfFile->wrf_temp()[lt] +
         coeff_interp[1]*wrfFile->wrf_temp()[rt];
   }

  amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, wrfInterptheta.begin(),
      wrfInterptheta.end(), m_wrf_theta_vals.begin());

  amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, tavg.line_average().begin(),
      tavg.line_average().end(), m_theta_vals.begin());

}

void ABLWrfForcingTemp::operator() (
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const
{

  const auto& dt = m_time.deltaT();

  const amrex::Real* wrftheta = m_wrf_theta_vals.data();
  const amrex::Real* thetavals = m_theta_vals.data();
  
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

    // // Compute Source term 
    src_term(i, j, k, 0) += (wrftheta[k] - thetavals[k]) / dt;

  });

}

}
}
}
