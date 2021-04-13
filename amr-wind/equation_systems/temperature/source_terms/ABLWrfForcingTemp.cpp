#include "amr-wind/equation_systems/temperature/source_terms/ABLWrfForcingTemp.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Print.H"

namespace amr_wind {
namespace pde {
namespace temperature {

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

ABLWrfForcingTemp::ABLWrfForcingTemp(const CFDSim& sim)
    : m_time(sim.time()), m_mesh(sim.mesh())
{
    const auto& abl = sim.physics_manager().get<amr_wind::ABL>();
    abl.register_wrf_temp_forcing(this);

    amrex::ParmParse pp(identifier());
    pp.query("forcing_scheme", m_forcing_scheme);
    
    mean_temperature_init(
        abl.abl_statistics().theta_profile(), abl.abl_wrf_file());
}

ABLWrfForcingTemp::~ABLWrfForcingTemp() = default;

void ABLWrfForcingTemp::mean_temperature_init(
    const FieldPlaneAveraging& tavg, const ABLWRFfile& wrfFile)
{
    m_axis = tavg.axis();
    // The implementation depends the assumption that the ABL statistics class
    // computes statistics at the cell-centeres only on level 0. If this
    // assumption changes in future, the implementation will break... so put in
    // a check here to catch this.
    AMREX_ALWAYS_ASSERT(
        m_mesh.Geom(0).Domain().length(m_axis) ==
        static_cast<int>(tavg.line_centroids().size()));

    m_theta_ht.resize(tavg.line_centroids().size());
    m_theta_vals.resize(tavg.ncell_line());

    m_nht = tavg.line_centroids().size();
    m_zht.resize(m_nht);

    m_error_wrf_avg_theta.resize(tavg.ncell_line());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, tavg.line_centroids().begin(),
        tavg.line_centroids().end(), m_theta_ht.begin());

    std::copy(tavg.line_centroids().begin(), tavg.line_centroids().end(),
            m_zht.begin());

    m_wrf_theta_vals.resize(wrfFile.nheights());
    m_wrf_ht.resize(wrfFile.nheights());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, wrfFile.wrf_heights().begin(),
        wrfFile.wrf_heights().end(), m_wrf_ht.begin());

    if (amrex::toLower(m_forcing_scheme) == "indirect") {
      indirectForcingInit();
    }

}

void ABLWrfForcingTemp::indirectForcingInit()
{

  amrex::Print() <<  "In indirect" << "\n";
  amrex::Real scaleFact = 1e-3; 

  amrex::Array2D<amrex::Real,0,3,0,3> zTz;

  // Generate the matrix Z^T W Z
  for (int irow=0; irow < 4; irow++){
    for (int icol=0; icol< 4; icol++){
      
      zTz(irow, icol) = 0.0;

      for(int iht=0; iht< m_nht; iht++){
        zTz(irow,icol) = zTz(irow,icol) + std::pow(m_zht[iht]*scaleFact, (icol+irow)); 
      }
    }
  }
  // Invert the matrix Z^T W Z
  invertMat(zTz, m_im_zTz);

}

void ABLWrfForcingTemp::invertMat(const amrex::Array2D<amrex::Real,0,3,0,3>& m, amrex::Array2D<amrex::Real,0,3,0,3>& im)
{

  amrex::Real A2323 = m(2, 2) * m(3, 3) - m(2, 3) * m(3, 2);
  amrex::Real A1323 = m(2, 1) * m(3, 3) - m(2, 3) * m(3, 1);
  amrex::Real A1223 = m(2, 1) * m(3, 2) - m(2, 2) * m(3, 1);
  amrex::Real A0323 = m(2, 0) * m(3, 3) - m(2, 3) * m(3, 0);
  amrex::Real A0223 = m(2, 0) * m(3, 2) - m(2, 2) * m(3, 0);
  amrex::Real A0123 = m(2, 0) * m(3, 1) - m(2, 1) * m(3, 0);
  amrex::Real A2313 = m(1, 2) * m(3, 3) - m(1, 3) * m(3, 2);
  amrex::Real A1313 = m(1, 1) * m(3, 3) - m(1, 3) * m(3, 1);
  amrex::Real A1213 = m(1, 1) * m(3, 2) - m(1, 2) * m(3, 1);
  amrex::Real A2312 = m(1, 2) * m(2, 3) - m(1, 3) * m(2, 2);
  amrex::Real A1312 = m(1, 1) * m(2, 3) - m(1, 3) * m(2, 1);
  amrex::Real A1212 = m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1);
  amrex::Real A0313 = m(1, 0) * m(3, 3) - m(1, 3) * m(3, 0);
  amrex::Real A0213 = m(1, 0) * m(3, 2) - m(1, 2) * m(3, 0);
  amrex::Real A0312 = m(1, 0) * m(2, 3) - m(1, 3) * m(2, 0);
  amrex::Real A0212 = m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0);
  amrex::Real A0113 = m(1, 0) * m(3, 1) - m(1, 1) * m(3, 0);
  amrex::Real A0112 = m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0);

  amrex::Real det = m(0, 0) * ( m(1, 1) * A2323 - m(1, 2) * A1323 + m(1, 3) * A1223 )
      - m(0, 1) * ( m(1, 0) * A2323 - m(1, 2) * A0323 + m(1, 3) * A0223 )
      + m(0, 2) * ( m(1, 0) * A1323 - m(1, 1) * A0323 + m(1, 3) * A0123 )
      - m(0, 3) * ( m(1, 0) * A1223 - m(1, 1) * A0223 + m(1, 2) * A0123 );
  det = 1.0 / det;

  im(0, 0) = det *   ( m(1, 1) * A2323 - m(1, 2) * A1323 + m(1, 3) * A1223 );
  im(0, 1) = det * - ( m(0, 1) * A2323 - m(0, 2) * A1323 + m(0, 3) * A1223 );
  im(0, 2) = det *   ( m(0, 1) * A2313 - m(0, 2) * A1313 + m(0, 3) * A1213 );
  im(0, 3) = det * - ( m(0, 1) * A2312 - m(0, 2) * A1312 + m(0, 3) * A1212 );
  im(1, 0) = det * - ( m(1, 0) * A2323 - m(1, 2) * A0323 + m(1, 3) * A0223 );
  im(1, 1) = det *   ( m(0, 0) * A2323 - m(0, 2) * A0323 + m(0, 3) * A0223 );
  im(1, 2) = det * - ( m(0, 0) * A2313 - m(0, 2) * A0313 + m(0, 3) * A0213 );
  im(1, 3) = det *   ( m(0, 0) * A2312 - m(0, 2) * A0312 + m(0, 3) * A0212 );
  im(2, 0) = det *   ( m(1, 0) * A1323 - m(1, 1) * A0323 + m(1, 3) * A0123 );
  im(2, 1) = det * - ( m(0, 0) * A1323 - m(0, 1) * A0323 + m(0, 3) * A0123 );
  im(2, 2) = det *   ( m(0, 0) * A1313 - m(0, 1) * A0313 + m(0, 3) * A0113 );
  im(2, 3) = det * - ( m(0, 0) * A1312 - m(0, 1) * A0312 + m(0, 3) * A0112 );
  im(3, 0) = det * - ( m(1, 0) * A1223 - m(1, 1) * A0223 + m(1, 2) * A0123 );
  im(3, 1) = det *   ( m(0, 0) * A1223 - m(0, 1) * A0223 + m(0, 2) * A0123 );
  im(3, 2) = det * - ( m(0, 0) * A1213 - m(0, 1) * A0213 + m(0, 2) * A0113 );
  im(3, 3) = det *   ( m(0, 0) * A1212 - m(0, 1) * A0212 + m(0, 2) * A0112 );
}



amrex::Real ABLWrfForcingTemp::mean_temperature_heights(
    const FieldPlaneAveraging& tavg, std::unique_ptr<ABLWRFfile>& wrfFile)
{

    amrex::Real currtime;
    currtime = m_time.current_time();

    // First the index in time
    m_idx_time = closest_index(wrfFile->wrf_times(), currtime);

    amrex::Array<amrex::Real, 2> coeff_interp{{0.0, 0.0}};

    amrex::Real denom =
        wrfFile->wrf_times()[m_idx_time + 1] - wrfFile->wrf_times()[m_idx_time];

    coeff_interp[0] = (wrfFile->wrf_times()[m_idx_time + 1] - currtime) / denom;
    coeff_interp[1] = 1.0 - coeff_interp[0];

    amrex::Real interpTflux;

    interpTflux = coeff_interp[0] * wrfFile->wrf_tflux()[m_idx_time] +
                  coeff_interp[1] * wrfFile->wrf_tflux()[m_idx_time + 1];

    int num_wrf_ht = wrfFile->nheights();

    amrex::Vector<amrex::Real> wrfInterptheta(num_wrf_ht);

    for (int i = 0; i < num_wrf_ht; i++) {
        int lt = m_idx_time * num_wrf_ht + i;
        int rt = (m_idx_time + 1) * num_wrf_ht + i;

        wrfInterptheta[i] = coeff_interp[0] * wrfFile->wrf_temp()[lt] +
                            coeff_interp[1] * wrfFile->wrf_temp()[rt];
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, wrfInterptheta.begin(), wrfInterptheta.end(),
        m_wrf_theta_vals.begin());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, tavg.line_average().begin(),
        tavg.line_average().end(), m_theta_vals.begin());

    size_t n_levels = tavg.ncell_line();
    amrex::Vector<amrex::Real> error_T(n_levels);
    
    for (size_t i = 0 ; i < n_levels; i++) {
      error_T[i] = wrfInterptheta[i]-tavg.line_average()[i];
    }

  
 if (amrex::toLower(m_forcing_scheme) == "indirect") {
     amrex::Array<amrex::Real, 4> ezP_T;

     amrex::Real scaleFact = 1e-3; 

     for (int i = 0; i < 4; i++) {
         ezP_T[i] = 0.0;

         for (int ih = 0; ih < m_nht; ih++) {
             ezP_T[i] = ezP_T[i] + error_T[ih] * std::pow(m_zht[ih]*scaleFact, i);
         }
     }

     for (int i = 0; i < 4; i++) {
         m_poly_coeff_theta[i] = 0.0;
         for (int j = 0; j < 4; j++) {
             m_poly_coeff_theta[i] = m_poly_coeff_theta[i] + m_im_zTz(i, j) * ezP_T[j];
         }
     }

     for (size_t ih = 0; ih < n_levels; ih++) {
         error_T[ih] = 0.0;
         for (int j = 0; j < 4; j++) {
             error_T[ih] =
                 error_T[ih] + m_poly_coeff_theta[j] * std::pow(m_zht[ih]*scaleFact, j);
         }
      }
 }


  amrex::Gpu::copy(
     amrex::Gpu::hostToDevice, error_T.begin(), error_T.end(), m_error_wrf_avg_theta.begin());
 

 return interpTflux;
}

void ABLWrfForcingTemp::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{

    const auto& dt = m_time.deltaT();

    const amrex::Real* wrftheta = m_wrf_theta_vals.data();
    const amrex::Real* thetavals = m_theta_vals.data();
    const amrex::Real* theta_error_val = m_error_wrf_avg_theta.data();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        // // Compute Source term
        // src_term(i, j, k, 0) += (wrftheta[k] - thetavals[k]) / dt;

        src_term(i, j, k, 0) += (theta_error_val[k])*0.2 / dt;
    });
}

} // namespace temperature
} // namespace pde
} // namespace amr_wind
