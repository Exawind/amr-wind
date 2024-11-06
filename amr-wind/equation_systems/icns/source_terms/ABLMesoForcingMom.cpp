#include "amr-wind/equation_systems/icns/source_terms/ABLMesoForcingMom.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/utilities/index_operations.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Print.H"
#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>
#include <iomanip>

namespace amr_wind::pde::icns {

ABLMesoForcingMom::ABLMesoForcingMom(const CFDSim& sim)
    : ABLMesoscaleForcing(sim, identifier())
{
    const auto& abl = sim.physics_manager().get<amr_wind::ABL>();
    abl.register_meso_mom_forcing(this);
    abl.abl_statistics().register_meso_mom_forcing(this);

    if (!abl.abl_meso_file().is_tendency_forcing()) {
        m_tendency = false;
        mean_velocity_init(
            abl.abl_statistics().vel_profile(), abl.abl_meso_file());
    } else {
        m_tendency = true;
        mean_velocity_init(abl.abl_meso_file());
    }

    if ((amrex::toLower(m_forcing_scheme) == "indirect") &&
        !m_update_transition_height) {
        indirect_forcing_init(); // do this once
    }
}

ABLMesoForcingMom::~ABLMesoForcingMom() = default;

void ABLMesoForcingMom::mean_velocity_init(const ABLMesoscaleInput& ncfile)
{

    const int num_meso_ht = ncfile.nheights();
    m_error_meso_avg_U.resize(num_meso_ht);
    m_error_meso_avg_V.resize(num_meso_ht);
    m_meso_u_vals.resize(num_meso_ht);
    m_meso_v_vals.resize(num_meso_ht);
    m_meso_ht.resize(num_meso_ht);
    m_err_U.resize(num_meso_ht);
    m_err_V.resize(num_meso_ht);

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, ncfile.meso_heights().begin(),
        ncfile.meso_heights().end(), m_meso_ht.begin());
}

void ABLMesoForcingMom::mean_velocity_init(
    const VelPlaneAveragingFine& vavg, const ABLMesoscaleInput& ncfile)
{
    const int num_meso_ht = ncfile.nheights();
    m_nht = vavg.ncell_line();
    m_axis = vavg.axis();
    // The implementation depends the assumption that the ABL statistics class
    // computes statistics at the cell-centeres only on level 0. If this
    // assumption changes in future, the implementation will break... so put in
    // a check here to catch this.
    AMREX_ALWAYS_ASSERT(m_mesh.Geom(0).Domain().length(m_axis) == m_nht);

    m_zht.resize(m_nht);
    m_vavg_ht.resize(m_nht);
    m_error_meso_avg_U.resize(m_nht);
    m_error_meso_avg_V.resize(m_nht);
    m_err_U.resize(m_nht);
    m_err_V.resize(m_nht);

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, vavg.line_centroids().begin(),
        vavg.line_centroids().end(), m_vavg_ht.begin());

    std::copy(
        vavg.line_centroids().begin(), vavg.line_centroids().end(),
        m_zht.begin());

    m_meso_u_vals.resize(num_meso_ht);
    m_meso_v_vals.resize(num_meso_ht);
    m_meso_ht.resize(num_meso_ht);

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, ncfile.meso_heights().begin(),
        ncfile.meso_heights().end(), m_meso_ht.begin());
}

void ABLMesoForcingMom::mean_velocity_heights(
    std::unique_ptr<ABLMesoscaleInput> const& ncfile)
{
    if (m_forcing_scheme.empty()) {
        return;
    }

    amrex::Real currtime;
    currtime = m_time.current_time();

    // First the index in time
    m_idx_time = utils::closest_index_ubound(ncfile->meso_times(), currtime);

    amrex::Array<amrex::Real, 2> coeff_interp{0.0, 0.0};

    amrex::Real denom =
        ncfile->meso_times()[m_idx_time + 1] - ncfile->meso_times()[m_idx_time];

    coeff_interp[0] = (ncfile->meso_times()[m_idx_time + 1] - currtime) / denom;
    coeff_interp[1] = 1.0 - coeff_interp[0];

    int num_meso_ht = ncfile->nheights();

    amrex::Vector<amrex::Real> time_interpolated_u(num_meso_ht);
    amrex::Vector<amrex::Real> time_interpolated_v(num_meso_ht);

    for (int i = 0; i < num_meso_ht; i++) {
        const int lt = m_idx_time * num_meso_ht + i;
        const int rt = (m_idx_time + 1) * num_meso_ht + i;

        time_interpolated_u[i] = coeff_interp[0] * ncfile->meso_u()[lt] +
                                 coeff_interp[1] * ncfile->meso_u()[rt];

        time_interpolated_v[i] = coeff_interp[0] * ncfile->meso_v()[lt] +
                                 coeff_interp[1] * ncfile->meso_v()[rt];
    }

    for (int ih = 0; ih < num_meso_ht; ih++) {
        m_err_U[ih] = time_interpolated_u[ih];
        m_err_V[ih] = time_interpolated_v[ih];
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, time_interpolated_u.begin(),
        time_interpolated_u.end(), m_error_meso_avg_U.begin());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, time_interpolated_v.begin(),
        time_interpolated_v.end(), m_error_meso_avg_V.begin());
}

void ABLMesoForcingMom::mean_velocity_heights(
    const VelPlaneAveragingFine& vavg,
    std::unique_ptr<ABLMesoscaleInput> const& ncfile)
{
    if (m_forcing_scheme.empty()) {
        return;
    }

    amrex::Real currtime;
    currtime = m_time.current_time();
    const auto& dt = m_time.delta_t();

    // First the index in time
    m_idx_time = utils::closest_index_ubound(ncfile->meso_times(), currtime);

    amrex::Array<amrex::Real, 2> coeff_interp{0.0, 0.0};

    amrex::Real denom =
        ncfile->meso_times()[m_idx_time + 1] - ncfile->meso_times()[m_idx_time];

    coeff_interp[0] = (ncfile->meso_times()[m_idx_time + 1] - currtime) / denom;
    coeff_interp[1] = 1.0 - coeff_interp[0];

    const int num_meso_ht = ncfile->nheights();

    amrex::Vector<amrex::Real> time_interpolated_u(num_meso_ht);
    amrex::Vector<amrex::Real> time_interpolated_v(num_meso_ht);

    for (int i = 0; i < num_meso_ht; i++) {
        const int lt = m_idx_time * num_meso_ht + i;
        const int rt = (m_idx_time + 1) * num_meso_ht + i;

        time_interpolated_u[i] = coeff_interp[0] * ncfile->meso_u()[lt] +
                                 coeff_interp[1] * ncfile->meso_u()[rt];

        time_interpolated_v[i] = coeff_interp[0] * ncfile->meso_v()[lt] +
                                 coeff_interp[1] * ncfile->meso_v()[rt];
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, time_interpolated_u.begin(),
        time_interpolated_u.end(), m_meso_u_vals.begin());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, time_interpolated_v.begin(),
        time_interpolated_v.end(), m_meso_v_vals.begin());

    const int numcomp = vavg.ncomp();

    amrex::Vector<amrex::Real> error_U(m_nht);
    amrex::Vector<amrex::Real> error_V(m_nht);

    const auto& vavg_lc = vavg.line_centroids();
    amrex::Vector<amrex::Real> meso_ht(num_meso_ht);
    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, m_meso_ht.begin(), m_meso_ht.end(),
        meso_ht.begin());
    for (int i = 0; i < m_nht; i++) {
        const amrex::Real height_interpolated_u =
            amr_wind::interp::linear(meso_ht, time_interpolated_u, vavg_lc[i]);
        const amrex::Real height_interpolated_v =
            amr_wind::interp::linear(meso_ht, time_interpolated_v, vavg_lc[i]);
        error_U[i] = height_interpolated_u -
                     vavg.line_average()[static_cast<int>(numcomp * i)];
        error_V[i] = height_interpolated_v -
                     vavg.line_average()[static_cast<int>(numcomp * i + 1)];
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

        amrex::Array<amrex::Real, 4> ezP_U;
        amrex::Array<amrex::Real, 4> ezP_V;

        // form Z^T W y
        for (int i = 0; i < 4; i++) {
            ezP_U[i] = 0.0;
            ezP_V[i] = 0.0;

            for (int ih = 0; ih < m_nht; ih++) {
                ezP_U[i] = ezP_U[i] + error_U[ih] * m_W[i] *
                                          std::pow(m_zht[ih] * m_scaleFact, i);
                ezP_V[i] = ezP_V[i] + error_V[ih] * m_W[i] *
                                          std::pow(m_zht[ih] * m_scaleFact, i);
            }
        }

        for (int i = 0; i < 4; i++) {
            m_poly_coeff_U[i] = 0.0;
            m_poly_coeff_V[i] = 0.0;
            for (int j = 0; j < 4; j++) {
                m_poly_coeff_U[i] =
                    m_poly_coeff_U[i] + m_im_zTz(i, j) * ezP_U[j];
                m_poly_coeff_V[i] =
                    m_poly_coeff_V[i] + m_im_zTz(i, j) * ezP_V[j];
            }
        }

        if (m_debug) {
            amrex::Print() << "direct vs indirect velocity error profile"
                           << std::endl;
        }
        amrex::Vector<amrex::Real> error_U_direct(m_nht);
        amrex::Vector<amrex::Real> error_V_direct(m_nht);
        for (int ih = 0; ih < m_nht; ih++) {
            error_U_direct[ih] = error_U[ih];
            error_V_direct[ih] = error_V[ih];
            error_U[ih] = 0.0;
            error_V[ih] = 0.0;
            for (int j = 0; j < 4; j++) {
                error_U[ih] =
                    error_U[ih] +
                    m_poly_coeff_U[j] * std::pow(m_zht[ih] * m_scaleFact, j);
                error_V[ih] =
                    error_V[ih] +
                    m_poly_coeff_V[j] * std::pow(m_zht[ih] * m_scaleFact, j);
            }

            if (m_debug) {
                amrex::Print() << m_zht[ih] << " " << error_U_direct[ih] << " "
                               << error_U[ih] << " " << error_V_direct[ih]
                               << " " << error_V[ih] << std::endl;
            }
        }

        if (amrex::toLower(m_forcing_transition) == "indirecttodirect") {
            blend_forcings(error_U, error_U_direct, error_U);
            blend_forcings(error_V, error_V_direct, error_V);

            if (m_debug) {
                for (int ih = 0; ih < m_nht; ih++) {
                    amrex::Print() << m_zht[ih] << " " << error_U[ih] << " "
                                   << error_V[ih] << std::endl;
                }
            }
        }
    }

    if (forcing_to_constant()) {
        constant_forcing_transition(error_U);
        constant_forcing_transition(error_V);

        if (m_debug) {
            for (int ih = 0; ih < m_nht; ih++) {
                amrex::Print() << m_zht[ih] << " " << error_U[ih] << " "
                               << error_V[ih] << std::endl;
            }
        }
    }

    for (int ih = 0; ih < m_nht; ih++) {
        error_U[ih] = error_U[ih] * m_gain_coeff / dt;
        error_V[ih] = error_V[ih] * m_gain_coeff / dt;
        m_err_U[ih] = error_U[ih];
        m_err_V[ih] = error_V[ih];
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, error_U.begin(), error_U.end(),
        m_error_meso_avg_U.begin());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, error_V.begin(), error_V.end(),
        m_error_meso_avg_V.begin());
}

void ABLMesoForcingMom::operator()(
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

    // The z values corresponding to the forcing values is either the number of
    // points in the netcdf input (for tendency) or the size of the plane
    // averaged velocities (non tendency)
    const amrex::Real* vheights_begin =
        (m_tendency) ? m_meso_ht.data() : m_vavg_ht.data();
    const amrex::Real* vheights_end = (m_tendency)
                                          ? m_meso_ht.data() + m_meso_ht.size()
                                          : m_vavg_ht.data() + m_vavg_ht.size();
    const amrex::Real* u_error_val = m_error_meso_avg_U.data();
    const amrex::Real* v_error_val = m_error_meso_avg_V.data();
    const int idir = (int)m_axis;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::IntVect iv(i, j, k);
        const amrex::Real ht = problo[idir] + (iv[idir] + 0.5) * dx[idir];
        const amrex::Real u_err = amr_wind::interp::linear(
            vheights_begin, vheights_end, u_error_val, ht);
        const amrex::Real v_err = amr_wind::interp::linear(
            vheights_begin, vheights_end, v_error_val, ht);

        // Compute Source term
        src_term(i, j, k, 0) += u_err;
        src_term(i, j, k, 1) += v_err;

        // No forcing in z-direction
    });
}

} // namespace amr_wind::pde::icns
