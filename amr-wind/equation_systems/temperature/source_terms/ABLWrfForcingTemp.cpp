#include "amr-wind/equation_systems/temperature/source_terms/ABLWrfForcingTemp.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Print.H"
#include "amr-wind/wind_energy/ABLWrf.H"
#include <AMReX_REAL.H>
#include <memory>

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
    : ABLWrfForcing(sim,identifier())
{
    const auto& abl = sim.physics_manager().get<amr_wind::ABL>();
    abl.register_wrf_temp_forcing(this);
    abl.abl_statistics().register_wrf_forcing_temp(this);

    if (!abl.abl_wrf_file().is_wrf_tendency_forcing()) {
      mean_temperature_init(
          abl.abl_statistics().theta_profile(), abl.abl_wrf_file());
    } else{
      mean_temperature_init(abl.abl_wrf_file());
    }

    if ((amrex::toLower(m_forcing_scheme) == "indirect") && !m_update_transition_height)
        indirectForcingInit(); // do this once
}

ABLWrfForcingTemp::~ABLWrfForcingTemp() = default;
void ABLWrfForcingTemp::mean_temperature_init(const ABLWRFfile& wrfFile){

  m_error_wrf_avg_theta.resize(wrfFile.nheights());
  m_wrf_theta_vals.resize(wrfFile.nheights());
  m_wrf_ht.resize(wrfFile.nheights());
  m_err_Theta.resize(wrfFile.nheights());

  amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, wrfFile.wrf_heights().begin(),
      wrfFile.wrf_heights().end(), m_wrf_ht.begin());

}
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

    m_err_Theta.resize(tavg.ncell_line());

    m_error_wrf_avg_theta.resize(tavg.ncell_line());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, tavg.line_centroids().begin(),
        tavg.line_centroids().end(), m_theta_ht.begin());

    std::copy(
        tavg.line_centroids().begin(), tavg.line_centroids().end(),
        m_zht.begin());

    m_wrf_theta_vals.resize(wrfFile.nheights());
    m_wrf_ht.resize(wrfFile.nheights());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, wrfFile.wrf_heights().begin(),
        wrfFile.wrf_heights().end(), m_wrf_ht.begin());
}

amrex::Real ABLWrfForcingTemp::mean_temperature_heights(std::unique_ptr<ABLWRFfile>& wrfFile){

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
        amrex::Gpu::hostToDevice, wrfInterptheta.begin(), wrfInterptheta.end(),m_error_wrf_avg_theta.begin());

    for (int ih = 0; ih < num_wrf_ht; ih++) {
      m_err_Theta[ih] = wrfInterptheta[ih];
    }

    return interpTflux;
  
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

    for (size_t i = 0; i < n_levels; i++) {
        error_T[i] = wrfInterptheta[i] - tavg.line_average()[i];
    }

    if (amrex::toLower(m_forcing_scheme) == "indirect") {
        if (m_update_transition_height) {
            //m_transition_height = 
            setTransitionWeighting();
            indirectForcingInit();
        }

        amrex::Array<amrex::Real, 4> ezP_T;

        // form Z^T W y
        for (int i = 0; i < 4; i++) {
            ezP_T[i] = 0.0;

            for (int ih = 0; ih < m_nht; ih++) {
                ezP_T[i] =
                    ezP_T[i] + error_T[ih] * m_W[ih] * std::pow(m_zht[ih] * m_scaleFact, i);
            }
        }

        for (int i = 0; i < 4; i++) {
            m_poly_coeff_theta[i] = 0.0;
            for (int j = 0; j < 4; j++) {
                m_poly_coeff_theta[i] =
                    m_poly_coeff_theta[i] + m_im_zTz(i, j) * ezP_T[j];
            }
        }

        if (m_debug)
            amrex::Print() << "direct vs indirect temperature error profile" << std::endl;
        amrex::Vector<amrex::Real> error_T_direct(n_levels);
        for (size_t ih = 0; ih < n_levels; ih++) {
            error_T_direct[ih] = error_T[ih];
            error_T[ih] = 0.0;
            for (int j = 0; j < 4; j++) {
                error_T[ih] =
                    error_T[ih] +
                    m_poly_coeff_theta[j] * std::pow(m_zht[ih] * m_scaleFact, j);
            }

            if (m_debug)
                amrex::Print() << m_zht[ih] << " " << error_T_direct[ih] << " " << error_T[ih] << std::endl;
        }

        if (amrex::toLower(m_forcing_transition) == "indirecttodirect") {
            blendForcings(error_T, error_T_direct, error_T);

            if (m_debug) {
                for (size_t ih=0; ih < n_levels; ih++) {
                    amrex::Print() << m_zht[ih] << " " << error_T[ih] << std::endl;
                }
            }
        }
    }

    if (forcingToConstant()) {
        constantForcingTransition(error_T);

        if (m_debug) {
            for (size_t ih=0; ih < n_levels; ih++) {
                amrex::Print() << m_zht[ih] << " " << error_T[ih] << std::endl;
            }
        }
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, error_T.begin(), error_T.end(),
        m_error_wrf_avg_theta.begin());

    for (size_t ih = 0; ih < n_levels; ih++) {
      m_err_Theta[ih] = error_T[ih]*m_gain_coeff;
    }

    return interpTflux;
}

void ABLWrfForcingTemp::operator()(
    const int lev,
    const amrex::MFIter&,
    const amrex::Box& bx,
    const FieldState,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& dt = m_time.deltaT();
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();

    const int idir = m_axis;
    const int nh_max =  m_wrf_ht.size() - 2; 
    const int lp1 = lev + 1;
    const amrex::Real* theights =  m_wrf_ht.data();
    const amrex::Real* theta_error_val = m_error_wrf_avg_theta.data();
    const amrex::Real kcoeff = m_gain_coeff;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::IntVect iv(i, j, k);
        const amrex::Real ht = problo[idir] + (iv[idir] + 0.5) * dx[idir];
        const int il = amrex::min(k / lp1, nh_max);
        const int ir = il + 1;
        amrex::Real temp;

        temp =
            theta_error_val[il] + ((theta_error_val[ir] - theta_error_val[il]) /
                                   (theights[ir] - theights[il])) *
                                      (ht - theights[il]);

        // Compute Source term
        // src_term(i, j, k, 0) += (theta_error_val[k]) * kcoeff / dt;
        src_term(i, j, k, 0) += temp * kcoeff / dt;
    });
}

} // namespace temperature
} // namespace pde
} // namespace amr_wind
