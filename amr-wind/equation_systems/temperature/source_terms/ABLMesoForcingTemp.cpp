#include "amr-wind/equation_systems/temperature/source_terms/ABLMesoForcingTemp.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/index_operations.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Print.H"
#include <AMReX_REAL.H>
#include <memory>

namespace amr_wind::pde::temperature {

ABLMesoForcingTemp::ABLMesoForcingTemp(const CFDSim& sim)
    : ABLMesoscaleForcing(sim, identifier())
{
    const auto& abl = sim.physics_manager().get<amr_wind::ABL>();
    abl.register_meso_temp_forcing(this);
    abl.abl_statistics().register_meso_temp_forcing(this);

    if (!abl.abl_meso_file().is_tendency_forcing()) {
        mean_temperature_init(
            abl.abl_statistics().theta_profile_fine(), abl.abl_meso_file());
    } else {
        mean_temperature_init(abl.abl_meso_file());
    }

    if ((amrex::toLower(m_forcing_scheme) == "indirect") &&
        !m_update_transition_height) {
        indirectForcingInit(); // do this once
    }
}

ABLMesoForcingTemp::~ABLMesoForcingTemp() = default;

void ABLMesoForcingTemp::mean_temperature_init(const ABLMesoscaleInput& ncfile)
{

    m_error_meso_avg_theta.resize(ncfile.nheights());
    m_meso_theta_vals.resize(ncfile.nheights());
    m_meso_ht.resize(ncfile.nheights());
    m_err_Theta.resize(ncfile.nheights());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, ncfile.meso_heights().begin(),
        ncfile.meso_heights().end(), m_meso_ht.begin());
}
void ABLMesoForcingTemp::mean_temperature_init(
    const FieldPlaneAveragingFine& tavg, const ABLMesoscaleInput& ncfile)
{
    m_axis = tavg.axis();
    // The implementation depends the assumption that the ABL statistics class
    // computes statistics at the cell-centeres only on level 0. If this
    // assumption changes in future, the implementation will break... so put in
    // a check here to catch this.
    AMREX_ALWAYS_ASSERT(
        m_mesh.Geom(0).Domain().length(m_axis) ==
        static_cast<int>(tavg.line_centroids().size()));

    m_nht = static_cast<int>(tavg.line_centroids().size());
    m_zht.resize(m_nht);

    m_theta_ht.resize(tavg.line_centroids().size());
    m_theta_vals.resize(tavg.ncell_line());

    m_err_Theta.resize(tavg.ncell_line());

    m_error_meso_avg_theta.resize(tavg.ncell_line());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, tavg.line_centroids().begin(),
        tavg.line_centroids().end(), m_theta_ht.begin());

    std::copy(
        tavg.line_centroids().begin(), tavg.line_centroids().end(),
        m_zht.begin());

    m_meso_theta_vals.resize(ncfile.nheights());
    m_meso_ht.resize(ncfile.nheights());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, ncfile.meso_heights().begin(),
        ncfile.meso_heights().end(), m_meso_ht.begin());
}

amrex::Real ABLMesoForcingTemp::mean_temperature_heights(
    std::unique_ptr<ABLMesoscaleInput> const& ncfile)
{

    amrex::Real currtime;
    currtime = m_time.current_time();

    // First the index in time
    m_idx_time = closest_index(ncfile->meso_times(), currtime);

    amrex::Array<amrex::Real, 2> coeff_interp{0.0, 0.0};

    amrex::Real denom =
        ncfile->meso_times()[m_idx_time + 1] - ncfile->meso_times()[m_idx_time];

    coeff_interp[0] = (ncfile->meso_times()[m_idx_time + 1] - currtime) / denom;
    coeff_interp[1] = 1.0 - coeff_interp[0];

    amrex::Real interpTflux;

    interpTflux = coeff_interp[0] * ncfile->meso_tflux()[m_idx_time] +
                  coeff_interp[1] * ncfile->meso_tflux()[m_idx_time + 1];

    if (m_forcing_scheme.empty()) {
        // no temperature profile assimilation
        return interpTflux;
    }

    int num_meso_ht = ncfile->nheights();

    amrex::Vector<amrex::Real> mesoInterptheta(num_meso_ht);

    for (int i = 0; i < num_meso_ht; i++) {
        int lt = m_idx_time * num_meso_ht + i;
        int rt = (m_idx_time + 1) * num_meso_ht + i;

        mesoInterptheta[i] = coeff_interp[0] * ncfile->meso_temp()[lt] +
                             coeff_interp[1] * ncfile->meso_temp()[rt];
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, mesoInterptheta.begin(),
        mesoInterptheta.end(), m_error_meso_avg_theta.begin());

    for (int ih = 0; ih < num_meso_ht; ih++) {
        m_err_Theta[ih] = mesoInterptheta[ih];
    }

    return interpTflux;
}

amrex::Real ABLMesoForcingTemp::mean_temperature_heights(
    const FieldPlaneAveragingFine& tavg,
    std::unique_ptr<ABLMesoscaleInput> const& ncfile)
{

    amrex::Real currtime;
    currtime = m_time.current_time();

    // First the index in time
    m_idx_time = closest_index(ncfile->meso_times(), currtime);

    amrex::Array<amrex::Real, 2> coeff_interp{0.0, 0.0};

    amrex::Real denom =
        ncfile->meso_times()[m_idx_time + 1] - ncfile->meso_times()[m_idx_time];

    coeff_interp[0] = (ncfile->meso_times()[m_idx_time + 1] - currtime) / denom;
    coeff_interp[1] = 1.0 - coeff_interp[0];

    amrex::Real interpTflux;

    interpTflux = coeff_interp[0] * ncfile->meso_tflux()[m_idx_time] +
                  coeff_interp[1] * ncfile->meso_tflux()[m_idx_time + 1];

    if (m_forcing_scheme.empty()) {
        // no temperature profile assimilation
        return interpTflux;
    }

    int num_meso_ht = ncfile->nheights();

    amrex::Vector<amrex::Real> mesoInterptheta(num_meso_ht);

    for (int i = 0; i < num_meso_ht; i++) {
        int lt = m_idx_time * num_meso_ht + i;
        int rt = (m_idx_time + 1) * num_meso_ht + i;

        mesoInterptheta[i] = coeff_interp[0] * ncfile->meso_temp()[lt] +
                             coeff_interp[1] * ncfile->meso_temp()[rt];
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, mesoInterptheta.begin(),
        mesoInterptheta.end(), m_meso_theta_vals.begin());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, tavg.line_average().begin(),
        tavg.line_average().end(), m_theta_vals.begin());

    size_t n_levels = tavg.ncell_line();
    amrex::Vector<amrex::Real> error_T(n_levels);

    for (size_t i = 0; i < n_levels; i++) {
        error_T[i] = mesoInterptheta[i] - tavg.line_average()[i];
    }

    if (amrex::toLower(m_forcing_scheme) == "indirect") {
        if (m_update_transition_height) {
            // possible unexpected behaviors, as described in
            // ec5eb95c6ca853ce0fea8488e3f2515a2d6374e7
            m_transition_height =
                coeff_interp[0] * ncfile->meso_transition_height()[m_idx_time] +
                coeff_interp[1] *
                    ncfile->meso_transition_height()[m_idx_time + 1];
            amrex::Print() << "current transition height = "
                           << m_transition_height << std::endl;

            setTransitionWeighting();
            indirectForcingInit();
        }

        amrex::Array<amrex::Real, 4> ezP_T;

        // form Z^T W y
        for (int i = 0; i < 4; i++) {
            ezP_T[i] = 0.0;

            for (int ih = 0; ih < m_nht; ih++) {
                ezP_T[i] = ezP_T[i] + error_T[ih] * m_W[ih] *
                                          std::pow(m_zht[ih] * m_scaleFact, i);
            }
        }

        for (int i = 0; i < 4; i++) {
            m_poly_coeff_theta[i] = 0.0;
            for (int j = 0; j < 4; j++) {
                m_poly_coeff_theta[i] =
                    m_poly_coeff_theta[i] + m_im_zTz(i, j) * ezP_T[j];
            }
        }

        if (m_debug) {
            amrex::Print() << "direct vs indirect temperature error profile"
                           << std::endl;
        }
        amrex::Vector<amrex::Real> error_T_direct(n_levels);
        for (size_t ih = 0; ih < n_levels; ih++) {
            error_T_direct[ih] = error_T[ih];
            error_T[ih] = 0.0;
            for (int j = 0; j < 4; j++) {
                error_T[ih] =
                    error_T[ih] + m_poly_coeff_theta[j] *
                                      std::pow(m_zht[ih] * m_scaleFact, j);
            }

            if (m_debug) {
                amrex::Print() << m_zht[ih] << " " << error_T_direct[ih] << " "
                               << error_T[ih] << std::endl;
            }
        }

        if (amrex::toLower(m_forcing_transition) == "indirecttodirect") {
            blendForcings(error_T, error_T_direct, error_T);

            if (m_debug) {
                for (size_t ih = 0; ih < n_levels; ih++) {
                    amrex::Print()
                        << m_zht[ih] << " " << error_T[ih] << std::endl;
                }
            }
        }
    }

    if (forcingToConstant()) {
        constantForcingTransition(error_T);

        if (m_debug) {
            for (size_t ih = 0; ih < n_levels; ih++) {
                amrex::Print() << m_zht[ih] << " " << error_T[ih] << std::endl;
            }
        }
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, error_T.begin(), error_T.end(),
        m_error_meso_avg_theta.begin());

    for (size_t ih = 0; ih < n_levels; ih++) {
        m_err_Theta[ih] = error_T[ih] * m_gain_coeff;
    }

    return interpTflux;
}

void ABLMesoForcingTemp::operator()(
    const int lev,
    const amrex::MFIter& /*mfi*/,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{
    if (m_forcing_scheme.empty()) {
        return;
    }

    const auto& dt = m_time.deltaT();
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();

    const int nh_max = (int)m_meso_ht.size() - 2;
    const int lp1 = lev + 1;
    const amrex::Real* theights = m_meso_ht.data();
    const amrex::Real* theta_error_val = m_error_meso_avg_theta.data();
    const amrex::Real kcoeff = m_gain_coeff;
    const int idir = (int)m_axis;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::IntVect iv(i, j, k);
        const amrex::Real ht = problo[idir] + (iv[idir] + 0.5) * dx[idir];
        const int il = amrex::min(k / lp1, nh_max);
        const int ir = il + 1;
        amrex::Real theta_err =
            theta_error_val[il] + ((theta_error_val[ir] - theta_error_val[il]) /
                                   (theights[ir] - theights[il])) *
                                      (ht - theights[il]);

        // Compute Source term
        src_term(i, j, k, 0) += kcoeff * theta_err / dt;
    });
}

} // namespace amr_wind::pde::temperature
