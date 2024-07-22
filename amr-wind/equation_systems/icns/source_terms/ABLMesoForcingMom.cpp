#include "amr-wind/equation_systems/icns/source_terms/ABLMesoForcingMom.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/utilities/index_operations.H"
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
        mean_velocity_init(
            abl.abl_statistics().vel_profile(), abl.abl_meso_file());
    } else {
        mean_velocity_init(abl.abl_meso_file());
    }

    if ((amrex::toLower(m_forcing_scheme) == "indirect") &&
        !m_update_transition_height) {
        indirectForcingInit(); // do this once
    }
}

ABLMesoForcingMom::~ABLMesoForcingMom() = default;

void ABLMesoForcingMom::mean_velocity_init(const ABLMesoscaleInput& ncfile)
{

    m_error_meso_avg_U.resize(ncfile.nheights());
    m_error_meso_avg_V.resize(ncfile.nheights());
    m_meso_u_vals.resize(ncfile.nheights());
    m_meso_v_vals.resize(ncfile.nheights());
    m_meso_ht.resize(ncfile.nheights());
    m_err_U.resize(ncfile.nheights());
    m_err_V.resize(ncfile.nheights());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, ncfile.meso_heights().begin(),
        ncfile.meso_heights().end(), m_meso_ht.begin());
}

void ABLMesoForcingMom::mean_velocity_init(
    const VelPlaneAveragingFine& vavg, const ABLMesoscaleInput& ncfile)
{
    m_axis = vavg.axis();
    // The implementation depends the assumption that the ABL statistics class
    // computes statistics at the cell-centeres only on level 0. If this
    // assumption changes in future, the implementation will break... so put in
    // a check here to catch this.
    AMREX_ALWAYS_ASSERT(
        m_mesh.Geom(0).Domain().length(m_axis) ==
        static_cast<int>(vavg.line_centroids().size()));

    m_nht = static_cast<int>(vavg.line_centroids().size());
    m_zht.resize(m_nht);

    m_velAvg_ht.resize(vavg.line_centroids().size());
    m_uAvg_vals.resize(vavg.ncell_line());
    m_vAvg_vals.resize(vavg.ncell_line());

    m_meso_avg_error.resize(vavg.ncell_line());

    m_error_meso_avg_U.resize(vavg.ncell_line());
    m_error_meso_avg_V.resize(vavg.ncell_line());

    m_err_U.resize(vavg.ncell_line());
    m_err_V.resize(vavg.ncell_line());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, vavg.line_centroids().begin(),
        vavg.line_centroids().end(), m_velAvg_ht.begin());

    std::copy(
        vavg.line_centroids().begin(), vavg.line_centroids().end(),
        m_zht.begin());

    m_meso_u_vals.resize(ncfile.nheights());
    m_meso_v_vals.resize(ncfile.nheights());
    m_meso_ht.resize(ncfile.nheights());

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
    m_idx_time = utils::closest_index(ncfile->meso_times(), currtime);

    amrex::Array<amrex::Real, 2> coeff_interp{0.0, 0.0};

    amrex::Real denom =
        ncfile->meso_times()[m_idx_time + 1] - ncfile->meso_times()[m_idx_time];

    coeff_interp[0] = (ncfile->meso_times()[m_idx_time + 1] - currtime) / denom;
    coeff_interp[1] = 1.0 - coeff_interp[0];

    int num_meso_ht = ncfile->nheights();

    amrex::Vector<amrex::Real> mesoInterpU(num_meso_ht);
    amrex::Vector<amrex::Real> mesoInterpV(num_meso_ht);

    for (int i = 0; i < num_meso_ht; i++) {
        int lt = m_idx_time * num_meso_ht + i;
        int rt = (m_idx_time + 1) * num_meso_ht + i;

        mesoInterpU[i] = coeff_interp[0] * ncfile->meso_u()[lt] +
                         coeff_interp[1] * ncfile->meso_u()[rt];

        mesoInterpV[i] = coeff_interp[0] * ncfile->meso_v()[lt] +
                         coeff_interp[1] * ncfile->meso_v()[rt];
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, mesoInterpU.begin(), mesoInterpU.end(),
        m_error_meso_avg_U.begin());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, mesoInterpV.begin(), mesoInterpV.end(),
        m_error_meso_avg_V.begin());

    for (int ih = 0; ih < num_meso_ht; ih++) {
        m_err_U[ih] = mesoInterpU[ih];
        m_err_V[ih] = mesoInterpV[ih];
    }
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

    // First the index in time
    m_idx_time = utils::closest_index(ncfile->meso_times(), currtime);

    amrex::Array<amrex::Real, 2> coeff_interp{0.0, 0.0};

    amrex::Real denom =
        ncfile->meso_times()[m_idx_time + 1] - ncfile->meso_times()[m_idx_time];

    coeff_interp[0] = (ncfile->meso_times()[m_idx_time + 1] - currtime) / denom;
    coeff_interp[1] = 1.0 - coeff_interp[0];

    int num_meso_ht = ncfile->nheights();

    amrex::Vector<amrex::Real> mesoInterpU(num_meso_ht);
    amrex::Vector<amrex::Real> mesoInterpV(num_meso_ht);

    for (int i = 0; i < num_meso_ht; i++) {
        int lt = m_idx_time * num_meso_ht + i;
        int rt = (m_idx_time + 1) * num_meso_ht + i;

        mesoInterpU[i] = coeff_interp[0] * ncfile->meso_u()[lt] +
                         coeff_interp[1] * ncfile->meso_u()[rt];

        mesoInterpV[i] = coeff_interp[0] * ncfile->meso_v()[lt] +
                         coeff_interp[1] * ncfile->meso_v()[rt];
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, mesoInterpU.begin(), mesoInterpU.end(),
        m_meso_u_vals.begin());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, mesoInterpV.begin(), mesoInterpV.end(),
        m_meso_v_vals.begin());

    // copy the spatially averaged velocity to GPU
    int numcomp = vavg.ncomp();
    size_t n_levels = vavg.ncell_line();
    amrex::Vector<amrex::Real> uStats(n_levels);
    amrex::Vector<amrex::Real> vStats(n_levels);
    for (size_t i = 0; i < n_levels; i++) {
        uStats[i] = vavg.line_average()[numcomp * i];
        vStats[i] = vavg.line_average()[numcomp * i + 1];
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, uStats.begin(), uStats.end(),
        m_uAvg_vals.begin());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, vStats.begin(), vStats.end(),
        m_vAvg_vals.begin());

    amrex::Vector<amrex::Real> error_U(n_levels);
    amrex::Vector<amrex::Real> error_V(n_levels);

    for (size_t i = 0; i < n_levels; i++) {
        error_U[i] = mesoInterpU[i] - uStats[i];
        error_V[i] = mesoInterpV[i] - vStats[i];
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
        amrex::Vector<amrex::Real> error_U_direct(n_levels);
        amrex::Vector<amrex::Real> error_V_direct(n_levels);
        for (size_t ih = 0; ih < n_levels; ih++) {
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
            blendForcings(error_U, error_U_direct, error_U);
            blendForcings(error_V, error_V_direct, error_V);

            if (m_debug) {
                for (size_t ih = 0; ih < n_levels; ih++) {
                    amrex::Print() << m_zht[ih] << " " << error_U[ih] << " "
                                   << error_V[ih] << std::endl;
                }
            }
        }
    }

    if (forcingToConstant()) {
        constantForcingTransition(error_U);
        constantForcingTransition(error_V);

        if (m_debug) {
            for (size_t ih = 0; ih < n_levels; ih++) {
                amrex::Print() << m_zht[ih] << " " << error_U[ih] << " "
                               << error_V[ih] << std::endl;
            }
        }
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, error_U.begin(), error_U.end(),
        m_error_meso_avg_U.begin());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, error_V.begin(), error_V.end(),
        m_error_meso_avg_V.begin());

    for (size_t ih = 0; ih < n_levels; ih++) {
        m_err_U[ih] = error_U[ih] * m_gain_coeff;
        m_err_V[ih] = error_V[ih] * m_gain_coeff;
    }
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

    const auto& dt = m_time.deltaT();
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();

    const int nh_max = (int)m_velAvg_ht.size() - 2;
    const int lp1 = lev + 1;
    const amrex::Real* vheights = m_meso_ht.data();
    const amrex::Real* u_error_val = m_error_meso_avg_U.data();
    const amrex::Real* v_error_val = m_error_meso_avg_V.data();
    const amrex::Real kcoeff = m_gain_coeff;
    const int idir = (int)m_axis;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::IntVect iv(i, j, k);
        const amrex::Real ht = problo[idir] + (iv[idir] + 0.5) * dx[idir];
        const int il = amrex::min(k / lp1, nh_max);
        const int ir = il + 1;

        amrex::Real u_err =
            u_error_val[il] + ((u_error_val[ir] - u_error_val[il]) /
                               (vheights[ir] - vheights[il])) *
                                  (ht - vheights[il]);

        amrex::Real v_err =
            v_error_val[il] + ((v_error_val[ir] - v_error_val[il]) /
                               (vheights[ir] - vheights[il])) *
                                  (ht - vheights[il]);

        // Compute Source term
        src_term(i, j, k, 0) += kcoeff * u_err / dt;
        src_term(i, j, k, 1) += kcoeff * v_err / dt;

        // No forcing in z-direction
    });
}

} // namespace amr_wind::pde::icns
