#include "amr-wind/equation_systems/icns/source_terms/ABLWrfForcingMom.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Print.H"
#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>
#include <iomanip>

namespace amr_wind {
namespace pde {
namespace icns {

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

ABLWrfForcingMom::ABLWrfForcingMom(const CFDSim& sim)
    : m_time(sim.time()), m_mesh(sim.mesh())
{

    const auto& abl = sim.physics_manager().get<amr_wind::ABL>();
    abl.register_mean_wrf_forcing(this);

    amrex::ParmParse pp(identifier());
    pp.query("forcing_scheme", m_forcing_scheme);

    mean_velocity_init(abl.abl_statistics().vel_profile(), abl.abl_wrf_file());

    pp.query("control_gain", m_gain_coeff);
}

ABLWrfForcingMom::~ABLWrfForcingMom() = default;

void ABLWrfForcingMom::mean_velocity_init(
    const VelPlaneAveraging& vavg, const ABLWRFfile& wrfFile)
{

    m_axis = vavg.axis();
    // The implementation depends the assumption that the ABL statistics class
    // computes statistics at the cell-centeres only on level 0. If this
    // assumption changes in future, the implementation will break... so put in
    // a check here to catch this.
    AMREX_ALWAYS_ASSERT(
        m_mesh.Geom(0).Domain().length(m_axis) ==
        static_cast<int>(vavg.line_centroids().size()));

    m_nht = vavg.line_centroids().size();
    m_zht.resize(m_nht);

    m_velAvg_ht.resize(vavg.line_centroids().size());
    m_uAvg_vals.resize(vavg.ncell_line());
    m_vAvg_vals.resize(vavg.ncell_line());

    m_wrf_avg_error.resize(vavg.ncell_line());

    m_error_wrf_avg_U.resize(vavg.ncell_line());
    m_error_wrf_avg_V.resize(vavg.ncell_line());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, vavg.line_centroids().begin(),
        vavg.line_centroids().end(), m_velAvg_ht.begin());

    std::copy(
        vavg.line_centroids().begin(), vavg.line_centroids().end(),
        m_zht.begin());

    m_wrf_u_vals.resize(wrfFile.nheights());
    m_wrf_v_vals.resize(wrfFile.nheights());
    m_wrf_ht.resize(wrfFile.nheights());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, wrfFile.wrf_heights().begin(),
        wrfFile.wrf_heights().end(), m_wrf_ht.begin());

    if (amrex::toLower(m_forcing_scheme) == "indirect") {
        indirectForcingInit();
    }

    m_nlevels = vavg.ncell_line();
}
void ABLWrfForcingMom::indirectForcingInit()
{

    amrex::Real scaleFact = 1e-3;

    amrex::Array2D<amrex::Real, 0, 3, 0, 3> zTz;

    // Generate the matrix Z^T W Z
    for (int irow = 0; irow < 4; irow++) {
        for (int icol = 0; icol < 4; icol++) {

            zTz(irow, icol) = 0.0;

            for (int iht = 0; iht < m_nht; iht++) {
                zTz(irow, icol) =
                    zTz(irow, icol) +
                    std::pow(m_zht[iht] * scaleFact, (icol + irow));
            }
        }
    }
    // Invert the matrix Z^T W Z
    invertMat(zTz, m_im_zTz);
}

void ABLWrfForcingMom::invertMat(
    const amrex::Array2D<amrex::Real, 0, 3, 0, 3>& m,
    amrex::Array2D<amrex::Real, 0, 3, 0, 3>& im)
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

    amrex::Real det =
        m(0, 0) * (m(1, 1) * A2323 - m(1, 2) * A1323 + m(1, 3) * A1223) -
        m(0, 1) * (m(1, 0) * A2323 - m(1, 2) * A0323 + m(1, 3) * A0223) +
        m(0, 2) * (m(1, 0) * A1323 - m(1, 1) * A0323 + m(1, 3) * A0123) -
        m(0, 3) * (m(1, 0) * A1223 - m(1, 1) * A0223 + m(1, 2) * A0123);
    det = 1.0 / det;

    im(0, 0) = det * (m(1, 1) * A2323 - m(1, 2) * A1323 + m(1, 3) * A1223);
    im(0, 1) = det * -(m(0, 1) * A2323 - m(0, 2) * A1323 + m(0, 3) * A1223);
    im(0, 2) = det * (m(0, 1) * A2313 - m(0, 2) * A1313 + m(0, 3) * A1213);
    im(0, 3) = det * -(m(0, 1) * A2312 - m(0, 2) * A1312 + m(0, 3) * A1212);
    im(1, 0) = det * -(m(1, 0) * A2323 - m(1, 2) * A0323 + m(1, 3) * A0223);
    im(1, 1) = det * (m(0, 0) * A2323 - m(0, 2) * A0323 + m(0, 3) * A0223);
    im(1, 2) = det * -(m(0, 0) * A2313 - m(0, 2) * A0313 + m(0, 3) * A0213);
    im(1, 3) = det * (m(0, 0) * A2312 - m(0, 2) * A0312 + m(0, 3) * A0212);
    im(2, 0) = det * (m(1, 0) * A1323 - m(1, 1) * A0323 + m(1, 3) * A0123);
    im(2, 1) = det * -(m(0, 0) * A1323 - m(0, 1) * A0323 + m(0, 3) * A0123);
    im(2, 2) = det * (m(0, 0) * A1313 - m(0, 1) * A0313 + m(0, 3) * A0113);
    im(2, 3) = det * -(m(0, 0) * A1312 - m(0, 1) * A0312 + m(0, 3) * A0112);
    im(3, 0) = det * -(m(1, 0) * A1223 - m(1, 1) * A0223 + m(1, 2) * A0123);
    im(3, 1) = det * (m(0, 0) * A1223 - m(0, 1) * A0223 + m(0, 2) * A0123);
    im(3, 2) = det * -(m(0, 0) * A1213 - m(0, 1) * A0213 + m(0, 2) * A0113);
    im(3, 3) = det * (m(0, 0) * A1212 - m(0, 1) * A0212 + m(0, 2) * A0112);
}

void ABLWrfForcingMom::mean_velocity_heights(
    const VelPlaneAveraging& vavg, std::unique_ptr<ABLWRFfile>& wrfFile)
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

    int num_wrf_ht = wrfFile->nheights();

    amrex::Vector<amrex::Real> wrfInterpU(num_wrf_ht);
    amrex::Vector<amrex::Real> wrfInterpV(num_wrf_ht);

    for (int i = 0; i < num_wrf_ht; i++) {
        int lt = m_idx_time * num_wrf_ht + i;
        int rt = (m_idx_time + 1) * num_wrf_ht + i;

        wrfInterpU[i] = coeff_interp[0] * wrfFile->wrf_u()[lt] +
                        coeff_interp[1] * wrfFile->wrf_u()[rt];

        wrfInterpV[i] = coeff_interp[0] * wrfFile->wrf_v()[lt] +
                        coeff_interp[1] * wrfFile->wrf_v()[rt];
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, wrfInterpU.begin(), wrfInterpU.end(),
        m_wrf_u_vals.begin());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, wrfInterpV.begin(), wrfInterpV.end(),
        m_wrf_v_vals.begin());

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
        error_U[i] = wrfInterpU[i] - uStats[i];
        error_V[i] = wrfInterpV[i] - vStats[i];
    }

    if (amrex::toLower(m_forcing_scheme) == "indirect") {
        amrex::Array<amrex::Real, 4> ezP_U;
        amrex::Array<amrex::Real, 4> ezP_V;

        amrex::Real scaleFact = 1e-3;

        for (int i = 0; i < 4; i++) {
            ezP_U[i] = 0.0;
            ezP_V[i] = 0.0;

            for (int ih = 0; ih < m_nht; ih++) {
                ezP_U[i] =
                    ezP_U[i] + error_U[ih] * std::pow(m_zht[ih] * scaleFact, i);
                ezP_V[i] =
                    ezP_V[i] + error_V[ih] * std::pow(m_zht[ih] * scaleFact, i);
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

        for (size_t ih = 0; ih < n_levels; ih++) {
            error_U[ih] = 0.0;
            error_V[ih] = 0.0;
            for (int j = 0; j < 4; j++) {
                error_U[ih] =
                    error_U[ih] +
                    m_poly_coeff_U[j] * std::pow(m_zht[ih] * scaleFact, j);
                error_V[ih] =
                    error_V[ih] +
                    m_poly_coeff_V[j] * std::pow(m_zht[ih] * scaleFact, j);
            }
        }
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, error_U.begin(), error_U.end(),
        m_error_wrf_avg_U.begin());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, error_V.begin(), error_V.end(),
        m_error_wrf_avg_V.begin());
}

void ABLWrfForcingMom::mean_velocity_heights(
    const amrex::Vector<amrex::Real>& uavg,
    const amrex::Vector<amrex::Real>& vavg,
    std::unique_ptr<ABLWRFfile>& wrfFile)
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

    int num_wrf_ht = wrfFile->nheights();

    amrex::Vector<amrex::Real> wrfInterpU(num_wrf_ht);
    amrex::Vector<amrex::Real> wrfInterpV(num_wrf_ht);

    for (int i = 0; i < num_wrf_ht; i++) {
        int lt = m_idx_time * num_wrf_ht + i;
        int rt = (m_idx_time + 1) * num_wrf_ht + i;

        wrfInterpU[i] = coeff_interp[0] * wrfFile->wrf_u()[lt] +
                        coeff_interp[1] * wrfFile->wrf_u()[rt];

        wrfInterpV[i] = coeff_interp[0] * wrfFile->wrf_v()[lt] +
                        coeff_interp[1] * wrfFile->wrf_v()[rt];
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, wrfInterpU.begin(), wrfInterpU.end(),
        m_wrf_u_vals.begin());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, wrfInterpV.begin(), wrfInterpV.end(),
        m_wrf_v_vals.begin());

    // copy the spatially averaged velocity to GPU
    amrex::Vector<amrex::Real> uStats(m_nlevels);
    amrex::Vector<amrex::Real> vStats(m_nlevels);
    for (size_t i = 0; i < m_nlevels; i++) {
        uStats[i] = uavg[i];
        vStats[i] = vavg[i];
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, uStats.begin(), uStats.end(),
        m_uAvg_vals.begin());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, vStats.begin(), vStats.end(),
        m_vAvg_vals.begin());

    amrex::Vector<amrex::Real> error_U(m_nlevels);
    amrex::Vector<amrex::Real> error_V(m_nlevels);

    for (size_t i = 0; i < m_nlevels; i++) {
        error_U[i] = wrfInterpU[i] - uStats[i];
        error_V[i] = wrfInterpV[i] - vStats[i];
    }

    if (amrex::toLower(m_forcing_scheme) == "indirect") {
        amrex::Array<amrex::Real, 4> ezP_U;
        amrex::Array<amrex::Real, 4> ezP_V;

        amrex::Real scaleFact = 1e-3;

        for (int i = 0; i < 4; i++) {
            ezP_U[i] = 0.0;
            ezP_V[i] = 0.0;

            for (int ih = 0; ih < m_nht; ih++) {
                ezP_U[i] =
                    ezP_U[i] + error_U[ih] * std::pow(m_zht[ih] * scaleFact, i);
                ezP_V[i] =
                    ezP_V[i] + error_V[ih] * std::pow(m_zht[ih] * scaleFact, i);
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

        for (size_t ih = 0; ih < m_nlevels; ih++) {
            error_U[ih] = 0.0;
            error_V[ih] = 0.0;
            for (int j = 0; j < 4; j++) {
                error_U[ih] =
                    error_U[ih] +
                    m_poly_coeff_U[j] * std::pow(m_zht[ih] * scaleFact, j);
                error_V[ih] =
                    error_V[ih] +
                    m_poly_coeff_V[j] * std::pow(m_zht[ih] * scaleFact, j);
            }
        }
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, error_U.begin(), error_U.end(),
        m_error_wrf_avg_U.begin());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, error_V.begin(), error_V.end(),
        m_error_wrf_avg_V.begin());
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

    const int idir = m_axis;
    const int nh_max = m_velAvg_ht.size() - 2;
    const int lp1 = lev + 1;
    const amrex::Real* vheights = m_velAvg_ht.data();
    const amrex::Real* u_error_val = m_error_wrf_avg_U.data();
    const amrex::Real* v_error_val = m_error_wrf_avg_V.data();
    const amrex::Real kcoeff = m_gain_coeff;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::IntVect iv(i, j, k);
        const amrex::Real ht = problo[idir] + (iv[idir] + 0.5) * dx[idir];
        const int il = amrex::min(k / lp1, nh_max);
        const int ir = il + 1;
        amrex::Real Utemp;
        amrex::Real Vtemp;

        Utemp = u_error_val[il] + ((u_error_val[ir] - u_error_val[il]) /
                                   (vheights[ir] - vheights[il])) *
                                      (ht - vheights[il]);

        Vtemp = v_error_val[il] + ((v_error_val[ir] - v_error_val[il]) /
                                   (vheights[ir] - vheights[il])) *
                                      (ht - vheights[il]);

        // // Compute Source term
        // src_term(i, j, k, 0) += u_error_val[k] * kcoeff / dt;
        // src_term(i, j, k, 1) += v_error_val[k] * kcoeff / dt;

        src_term(i, j, k, 0) += Utemp * kcoeff / dt;
        src_term(i, j, k, 1) += Vtemp * kcoeff / dt;

        // No forcing in z-direction
    });
}

} // namespace icns
} // namespace pde
} // namespace amr_wind
