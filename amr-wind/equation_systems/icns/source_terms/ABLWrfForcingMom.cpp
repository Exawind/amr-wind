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

}

ABLWrfForcingMom::~ABLWrfForcingMom() = default;

void ABLWrfForcingMom::mean_velocity_init(const VelPlaneAveraging& vavg, std::shared_ptr<ABLWRFfile>& wrfFile)
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

    m_wrf_u_vals.resize(wrfFile->nheights());
    m_wrf_v_vals.resize(wrfFile->nheights());
    m_wrf_ht.resize(wrfFile->nheights());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, wrfFile->wrf_heights().begin(),
        wrfFile->wrf_heights().end(), m_wrf_ht.begin());
    
    mean_velocity_heights(vavg, wrfFile);
}

void ABLWrfForcingMom::mean_velocity_heights(const VelPlaneAveraging& vavg, std::shared_ptr<ABLWRFfile>& wrfFile)
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

  amrex::Print() << "size of avg array = " << vavg.line_average().size() << "\n";

  // copy the spatially averaged velocity to GPU 
  // size_t n_levels = vavg.ncell_line();
  // amrex::Vector<amrex::Real> l_vec(n_levels);

  // vavg.line_average(0, l_vec);
  // amrex::Gpu::copy(
  //     amrex::Gpu::hostToDevice, l_vec.begin(), l_vec.end(), m_uAvg_vals.begin());

  // vavg.line_average(1, l_vec);
  // amrex::Gpu::copy(
  //     amrex::Gpu::hostToDevice, l_vec.begin(), l_vec.end(), m_vAvg_vals.begin());
      
}

void ABLWrfForcingMom::operator()(
    const int lev,
    const amrex::MFIter&,
    const amrex::Box& bx,
    const FieldState,
    const amrex::Array4<amrex::Real>& src_term) const
{
  const auto& dt = m_time.deltaT();
  const auto& problo = m_mesh.Geom(lev).ProbLoArray();
  const auto& dx = m_mesh.Geom(lev).CellSizeArray();

  const amrex::Real* wrfu = m_wrf_u_vals.data();
  const amrex::Real* wrfv = m_wrf_v_vals.data();

  const amrex::Real* velheights = m_velAvg_ht.data();
  const amrex::Real* uvals = m_uAvg_vals.data();
  const amrex::Real* vvals = m_vAvg_vals.data();

  const int idir = m_axis;
  const int lp1 = lev + 1;
  
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

    amrex::IntVect iv(i, j, k);
    const amrex::Real ht = problo[idir] + (iv[idir] + 0.5) * dx[idir];

    // First interpolate WRF output to the cell center
    // int idx_ht;  
    // // idx_ht = closest_index(m_wrf_ht, ht);

    // const int nh_max_wrf = m_wrf_ht.size()-2; 

    // int ilh = amrex::min(idx_ht, nh_max_wrf);
    // int irh = ilh + 1;

    // amrex::Real wrfIntp_u = wrfu[ilh] +
    //     ((wrfu[irh] - wrfu[ilh]) /
    //     (m_wrf_ht[irh] - m_wrf_ht[ilh])) *
    //     (ht - m_wrf_ht[ilh]);

    // amrex::Real wrfIntp_v = wrfv[ilh] +
    //     ((wrfv[irh] - wrfv[ilh]) /
    //      (m_wrf_ht[irh] - m_wrf_ht[ilh])) *
    //     (ht - m_wrf_ht[ilh]);

    // Read in cell center spatially averaged u and v velocity
    const int nh_max = m_velAvg_ht.size() - 2;
    const int il = amrex::min(k / lp1, nh_max);
    const int ir = il + 1;
    amrex::Real meanU, meanV;
    
    meanU = uvals[il] +
        ((uvals[ir] - uvals[il]) / (velheights[ir] - velheights[il])) *
        (ht - velheights[il]);

    meanV = vvals[il] +
        ((vvals[ir] - vvals[il]) / (velheights[ir] - velheights[il])) *
        (ht - velheights[il]);

    // Compute Source term 
    src_term(i, j, k, 0) += (wrfu[k] - meanU) / dt;
    src_term(i, j, k, 1) += (wrfv[k] - meanV) / dt;
    
    // No forcing in z-direction
  });

  
}
                                                        
} // namespace icns
} // namespace pde
} // namespace amr_wind
