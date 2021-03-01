#include "amr-wind/equation_systems/icns/source_terms/ABLWrfForcingMom.H"
#include "amr-wind/utilities/PlaneAveraging.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABL.H"


#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"

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

  amrex::Vector<amrex::Real> tmpHeights(m_nheight);
  
  ncf.var("heights").get(tmpHeights.data());
  ncf.var("times").get(m_wrf_time.data());

  m_wrf_mom.resize(m_nheight*m_ntime*2);

  ncf.var("wrf_momentum").get(m_wrf_mom.data());

  ncutils::NCFile::close();

  amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, tmpHeights.begin(),
      tmpHeights.end(), m_wrf_height.begin());


}

void ABLWrfForcingMom::mean_velocity_init(const FieldPlaneAveraging& vavg)
{
  // The implementation depends the assumption that the ABL statistics class
  // computes statistics at the cell-centeres only on level 0. If this
  // assumption changes in future, the implementation will break... so put in
  // a check here to catch this.
    AMREX_ALWAYS_ASSERT(
        m_mesh.Geom(0).Domain().length(m_axis) ==
        static_cast<int>(tavg.line_centroids().size()));
    m_velAvg_ht.resize(vavg.line_centroids().size());
    m_uAvg_vals.resize(vavg.ncell_line());
    m_vAvg_vals.resize(vavg.ncell_line());
    
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, vavg.line_centroids().begin(),
        vavg.line_centroids().end(), m_velAvg_ht.begin());

    m_wrf_u_vals.resize(m_nheight);
    m_wrf_v_vals.resize(m_nheight);
    
    mean_velocity_heights(vavg);
}

void ABLWrfForcingMom::mean_velocity_heights(const FieldPlaneAveraging& vavg)
{

  amrex::Real currtime;
  currtime = m_time.current_time();
  
  // First the index in time
  m_idx_time = closest_index(m_wrf_time, currtime);

  amrex::Array<amrex::Real, 2> coeff_interp{{0.0, 0.0}};

  coeff_interp[0] =  (m_wrf_time[m_idx_time+1] - m_time) /
      (m_wrf_time[m_idx_time+1] - m_wrf_time[m_idx_time]);
  coeff_interp[1] = 1.0 - coeff_interp[0];

  amrex::Vector<amrex::Real> wrfInterpU(m_nheight);
  amrex::Vector<amrex::Real> wrfInterpV(m_nheight);

  for (i = 0 ; i < m_nheight; i++)
  {
    int lt = m_idx_time*m_nheight*2 + i*2;
    int rt = (m_idx_time+1)*m_nheight*2 + i*2;
    wrfInterpU[i] = coeff_interp[0]*m_wrf_mom[lt] +
        coeff_interp[1]*m_wrf_mom[rt];
    
    wrfInterpV[i] = coeff_interp[0]*m_wrf_mom[lt+1] +
        coeff_interp[1]*m_wrf_mom[rt+1];
  }

  amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, wrfInterpU.begin(),
      wrfInterpU.end(), m_wrf_u_vals.begin());

  amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, wrfInterpV.begin(),
      wrfInterpV.end(), m_wrf_v_vals.begin());


  // copy the spatially averaged velocity to GPU 
  size_t n_levels = vavg.ncell_line();
  amrex::Vector<amrex::Real> l_vec(n_levels);
  vavg.line_average(0, l_vec);

  amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, l_vec.begin(), l_vec.end(), m_uAvg_vals.begin());

  vavg.line_average(1, l_vec);

  amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, l_vec.begin(), l_vec.end(), m_vAvg_vals.begin());
      
}

void ABLForcing::operator()(
    const int,
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

  const int nh_max = m_velAvg_ht.size() - 2;
  const int lp1 = lev + 1;
  

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

    amrex::IntVect iv(i, j, k);
    const amrex::Real ht = problo[idir] + (iv[idir] + 0.5) * dx[idir];

     // First interpolate WRF output to the cell center
    int idx_ht;  
    idx_ht = closest_index(m_wrf_height, ht);

    int ilh = amrex::min(idx_ht, m_nheight-2);
    int irh = ilh + 1;

    amrex::Real wrfIntp_u = wrfu[ilh] +
        ((wrfu[irh] - wrfu[ilh]) /
         (m_wrf_height[irh] - m_wrf_height[ilh])) *
        (ht - m_wrf_height[ilh]);

    amrex::Real wrfIntp_v = wrfv[ilh] +
        ((wrfv[irh] - wrfv[ilh]) /
         (m_wrf_height[irh] - m_wrf_height[ilh])) *
        (ht - m_wrf_height[ilh]);

    // Read in cell center spatially averaged u and v velocity
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
    src_term(i, j, k, 0) += (wrfIntp_u - meanU) / dt;
    src_term(i, j, k, 1) += (wrfIntp_v - meanV) / dt;
    
    // No forcing in z-direction
  });

  
}
                                                        
} // namespace icns
} // namespace pde
} // namespace amr_wind
