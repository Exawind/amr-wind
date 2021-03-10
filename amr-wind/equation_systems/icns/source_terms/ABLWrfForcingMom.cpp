#include "amr-wind/equation_systems/icns/source_terms/ABLWrfForcingMom.H"
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
}

ABLWrfForcingMom::ABLWrfForcingMom(const CFDSim& sim)
    : m_time(sim.time()), m_mesh(sim.mesh())
{

  const auto& abl = sim.physics_manager().get<amr_wind::ABL>();
  abl.register_mean_wrf_forcing(this);

  mean_velocity_init(abl.abl_statistics().vel_profile(), abl.abl_wrf_file());
  
}

ABLWrfForcingMom::~ABLWrfForcingMom() = default;

void ABLWrfForcingMom::mean_velocity_init(const VelPlaneAveraging& vavg, const ABLWRFfile& wrfFile)
{

  m_axis = vavg.axis();
  // The implementation depends the assumption that the ABL statistics class
  // computes statistics at the cell-centeres only on level 0. If this
  // assumption changes in future, the implementation will break... so put in
  // a check here to catch this.
  AMREX_ALWAYS_ASSERT(
      m_mesh.Geom(0).Domain().length(m_axis) ==
      static_cast<int>(vavg.line_centroids().size()));
  
  m_velAvg_ht.resize(vavg.line_centroids().size());
  m_uAvg_vals.resize(vavg.ncell_line());
  m_vAvg_vals.resize(vavg.ncell_line());
    
  amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, vavg.line_centroids().begin(),
      vavg.line_centroids().end(), m_velAvg_ht.begin());

  m_wrf_u_vals.resize(wrfFile.nheights());
  m_wrf_v_vals.resize(wrfFile.nheights());
  m_wrf_ht.resize(wrfFile.nheights());

  amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, wrfFile.wrf_heights().begin(),
      wrfFile.wrf_heights().end(), m_wrf_ht.begin());
    
}

void ABLWrfForcingMom::mean_velocity_heights(const VelPlaneAveraging& vavg, std::unique_ptr<ABLWRFfile>& wrfFile)
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
  
  amrex::Vector<amrex::Real> wrfInterpU(num_wrf_ht);
  amrex::Vector<amrex::Real> wrfInterpV(num_wrf_ht);

   for (int i = 0 ; i < num_wrf_ht; i++) {
     int lt = m_idx_time*num_wrf_ht + i;
     int rt = (m_idx_time+1)*num_wrf_ht + i;
    
     wrfInterpU[i] = coeff_interp[0]*wrfFile->wrf_u()[lt] +
         coeff_interp[1]*wrfFile->wrf_u()[rt];
    
     wrfInterpV[i] = coeff_interp[0]*wrfFile->wrf_v()[lt] +
         coeff_interp[1]*wrfFile->wrf_v()[rt];
   }

  amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, wrfInterpU.begin(),
      wrfInterpU.end(), m_wrf_u_vals.begin());

  amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, wrfInterpV.begin(),
      wrfInterpV.end(), m_wrf_v_vals.begin());

  // copy the spatially averaged velocity to GPU
  int numcomp = vavg.ncomp();
  size_t n_levels = vavg.ncell_line();
  amrex::Vector<amrex::Real> uStats(n_levels);
  amrex::Vector<amrex::Real> vStats(n_levels);
 for (int i = 0 ; i < n_levels; i++) {
   uStats[i] = vavg.line_average()[numcomp*i];
   vStats[i] = vavg.line_average()[numcomp*i + 1];
 }

 amrex::Gpu::copy(
     amrex::Gpu::hostToDevice, uStats.begin(), uStats.end(), m_uAvg_vals.begin());

 amrex::Gpu::copy(
     amrex::Gpu::hostToDevice, vStats.begin(), vStats.end(), m_vAvg_vals.begin());
      
}

void ABLWrfForcingMom::operator()(
    const int lev,
    const amrex::MFIter&,
    const amrex::Box& bx,
    const FieldState,
    const amrex::Array4<amrex::Real>& src_term) const
{
  const auto& dt = m_time.deltaT();

  const amrex::Real* wrfu = m_wrf_u_vals.data();
  const amrex::Real* wrfv = m_wrf_v_vals.data();

  const amrex::Real* velheights = m_velAvg_ht.data();
  const amrex::Real* uvals = m_uAvg_vals.data();
  const amrex::Real* vvals = m_vAvg_vals.data();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

    // // Compute Source term 
    src_term(i, j, k, 0) += (wrfu[k] - uvals[k]) / dt;
    src_term(i, j, k, 1) += (wrfv[k] - vvals[k]) / dt;
    
    // No forcing in z-direction

  });

  
}
                                                        
} // namespace icns
} // namespace pde
} // namespace amr_wind
