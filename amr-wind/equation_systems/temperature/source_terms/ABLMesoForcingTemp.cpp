#include "amr-wind/equation_systems/temperature/source_terms/ABLMesoForcingTemp.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/index_operations.H"
#include "amr-wind/utilities/linear_interpolation.H"
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
        indirect_forcing_init(); // do this once
    }
}

ABLMesoForcingTemp::~ABLMesoForcingTemp() = default;

void ABLMesoForcingTemp::mean_temperature_init(const ABLMesoscaleInput& ncfile)
{

    const int num_meso_ht = ncfile.nheights();
    m_error_meso_avg_theta.resize(num_meso_ht);
    m_meso_theta_vals.resize(num_meso_ht);
    m_meso_ht.resize(num_meso_ht);
    m_err_Theta.resize(num_meso_ht);

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, ncfile.meso_heights().begin(),
        ncfile.meso_heights().end(), m_meso_ht.begin());
}
void ABLMesoForcingTemp::mean_temperature_init(
    const FieldPlaneAveragingFine& tavg, const ABLMesoscaleInput& ncfile)
{
    const int num_meso_ht = ncfile.nheights();
    m_nht = tavg.ncell_line();

    m_axis = tavg.axis();
    // The implementation depends the assumption that the ABL statistics class
    // computes statistics at the cell-centeres only on level 0. If this
    // assumption changes in future, the implementation will break... so put in
    // a check here to catch this.
    AMREX_ALWAYS_ASSERT(m_mesh.Geom(0).Domain().length(m_axis) == m_nht);

    m_zht.resize(m_nht);

    m_theta_ht.resize(m_nht);
    m_theta_vals.resize(m_nht);

    m_err_Theta.resize(m_nht);

    m_error_meso_avg_theta.resize(m_nht);

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, tavg.line_centroids().begin(),
        tavg.line_centroids().end(), m_theta_ht.begin());

    std::copy(
        tavg.line_centroids().begin(), tavg.line_centroids().end(),
        m_zht.begin());

    m_meso_theta_vals.resize(num_meso_ht);
    m_meso_ht.resize(num_meso_ht);

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
    m_idx_time = utils::closest_index(ncfile->meso_times(), currtime);

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

    const int num_meso_ht = ncfile->nheights();

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
    const auto& dt = m_time.delta_t();

    // First the index in time
    m_idx_time = utils::closest_index(ncfile->meso_times(), currtime);

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

    const int num_meso_ht = ncfile->nheights();

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

    amrex::Vector<amrex::Real> error_T(m_nht);

    for (size_t i = 0; i < m_nht; i++) {
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

            set_transition_weighting();
            indirect_forcing_init();
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
        amrex::Vector<amrex::Real> error_T_direct(m_nht);
        for (size_t ih = 0; ih < m_nht; ih++) {
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
            blend_forcings(error_T, error_T_direct, error_T);

            if (m_debug) {
                for (size_t ih = 0; ih < m_nht; ih++) {
                    amrex::Print()
                        << m_zht[ih] << " " << error_T[ih] << std::endl;
                }
            }
        }
    }

    if (forcing_to_constant()) {
        constant_forcing_transition(error_T);

        if (m_debug) {
            for (size_t ih = 0; ih < m_nht; ih++) {
                amrex::Print() << m_zht[ih] << " " << error_T[ih] << std::endl;
            }
        }
    }

    for (size_t ih = 0; ih < m_nht; ih++) {
        error_T[ih] = error_T[ih] * m_gain_coeff / dt;
        m_err_Theta[ih] = error_T[ih];
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, error_T.begin(), error_T.end(),
        m_error_meso_avg_theta.begin());

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

    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();

    const amrex::Real* theights_begin = m_meso_ht.data();
    const amrex::Real* theights_end = m_meso_ht.data() + m_meso_ht.size();
    const amrex::Real* theta_error_val = m_error_meso_avg_theta.data();
    const int idir = (int)m_axis;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::IntVect iv(i, j, k);
        const amrex::Real ht = problo[idir] + (iv[idir] + 0.5) * dx[idir];
        const amrex::Real theta_err = amr_wind::interp::linear(
            theights_begin, theights_end, theta_error_val, ht);

        // Compute Source term
        src_term(i, j, k, 0) += theta_err;
    });
}

} // namespace amr_wind::pde::temperature
