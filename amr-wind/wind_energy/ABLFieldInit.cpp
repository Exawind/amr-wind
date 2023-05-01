#include <cmath>

#include "amr-wind/wind_energy/ABLFieldInit.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

ABLFieldInit::ABLFieldInit()
{
    amrex::ParmParse pp_abl("ABL");

    // Temperature variation as a function of height
    pp_abl.getarr("temperature_heights", m_theta_heights);
    pp_abl.getarr("temperature_values", m_theta_values);

    AMREX_ALWAYS_ASSERT(m_theta_heights.size() == m_theta_values.size());
    int num_theta_values = static_cast<int>(m_theta_heights.size());

    pp_abl.query("perturb_velocity", m_perturb_vel);
    pp_abl.query("perturb_ref_height", m_ref_height);
    pp_abl.query("Uperiods", m_Uperiods);
    pp_abl.query("Vperiods", m_Vperiods);
    pp_abl.query("deltaU", m_deltaU);
    pp_abl.query("deltaV", m_deltaV);

    pp_abl.query("perturb_temperature", m_perturb_theta);
    pp_abl.query("random_gauss_mean", m_theta_gauss_mean);
    pp_abl.query("random_gauss_var", m_theta_gauss_var);
    pp_abl.query("cutoff_height", m_theta_cutoff_height);
    pp_abl.query("theta_amplitude", m_deltaT);

    pp_abl.query("init_tke", m_tke_init);
    pp_abl.query("init_tke_beare_profile", m_tke_init_profile);
    pp_abl.query("init_tke_beare_factor", m_tke_init_factor);
    pp_abl.query("init_tke_cutoff_height", m_tke_cutoff_height);

    pp_abl.query("linear_profile", m_linear_profile);

    pp_abl.query("top_velocity", m_top_vel);
    pp_abl.query("bottom_velocity", m_bottom_vel);

    // TODO: Modify this to accept velocity as a function of height
    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.get("density", m_rho);

    amrex::ParmParse pp_forcing("ABLForcing");
    pp_forcing.query("velocity_timetable", m_vel_timetable);
    if (!m_vel_timetable.empty()) {
        std::ifstream ifh(m_vel_timetable, std::ios::in);
        if (!ifh.good()) {
            amrex::Abort("Cannot find input file: " + m_vel_timetable);
        }
        amrex::Real m_vel_time;
        amrex::Real m_vel_ang;
        ifh >> m_vel_time >> m_vel_speed >> m_vel_ang;
        m_vel_dir = utils::radians(m_vel_ang);
    } else {
        pp_incflo.getarr("velocity", m_vel);
    }

    m_thht_d.resize(num_theta_values);
    m_thvv_d.resize(num_theta_values);

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_theta_heights.begin(),
        m_theta_heights.end(), m_thht_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_theta_values.begin(), m_theta_values.end(),
        m_thvv_d.begin());
}

void ABLFieldInit::operator()(
    const amrex::Box& vbx,
    const amrex::Geometry& geom,
    const amrex::Array4<amrex::Real>& velocity,
    const amrex::Array4<amrex::Real>& density,
    const amrex::Array4<amrex::Real>& temperature) const
{
    const amrex::Real pi = M_PI;
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    const auto& probhi = geom.ProbHiArray();

    const bool perturb_vel = m_perturb_vel;

    const bool linear_profile = m_linear_profile;

    const amrex::Real rho_init = m_rho;

    const amrex::Real umean =
        !m_vel_timetable.empty() ? m_vel_speed * std::cos(m_vel_dir) : m_vel[0];
    const amrex::Real vmean =
        !m_vel_timetable.empty() ? m_vel_speed * std::sin(m_vel_dir) : m_vel[1];
    const amrex::Real wmean = !m_vel_timetable.empty() ? 0.0 : m_vel[2];

    const amrex::Real top_u_vel = m_top_vel[0];
    const amrex::Real top_v_vel = m_top_vel[1];
    const amrex::Real top_w_vel = m_top_vel[2];

    const amrex::Real bottom_u_vel = m_bottom_vel[0];
    const amrex::Real bottom_v_vel = m_bottom_vel[1];
    const amrex::Real bottom_w_vel = m_bottom_vel[2];

    const amrex::Real aval = m_Uperiods * 2.0 * pi / (probhi[1] - problo[1]);
    const amrex::Real bval = m_Vperiods * 2.0 * pi / (probhi[0] - problo[0]);
    const amrex::Real ufac = m_deltaU * std::exp(0.5) / m_ref_height;
    const amrex::Real vfac = m_deltaV * std::exp(0.5) / m_ref_height;
    const amrex::Real ref_height = m_ref_height;

    const int ntvals = static_cast<int>(m_theta_heights.size());
    const amrex::Real* th = m_thht_d.data();
    const amrex::Real* tv = m_thvv_d.data();

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
        const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

        density(i, j, k) = rho_init;
        // Mean velocity field
        velocity(i, j, k, 0) = umean;
        velocity(i, j, k, 1) = vmean;
        velocity(i, j, k, 2) = wmean;

        amrex::Real theta = tv[0];
        for (int iz = 0; iz < ntvals - 1; ++iz) {
            if ((z > th[iz]) && (z <= th[iz + 1])) {
                const amrex::Real slope =
                    (tv[iz + 1] - tv[iz]) / (th[iz + 1] - th[iz]);
                theta = tv[iz] + (z - th[iz]) * slope;
            }
        }

        temperature(i, j, k, 0) += theta;

        if (linear_profile) {
            velocity(i, j, k, 0) =
                bottom_u_vel +
                z * (top_u_vel - bottom_u_vel) / (probhi[2] - problo[2]);
            velocity(i, j, k, 1) =
                bottom_v_vel +
                z * (top_v_vel - bottom_v_vel) / (probhi[2] - problo[2]);
            velocity(i, j, k, 2) =
                bottom_w_vel +
                z * (top_w_vel - bottom_w_vel) / (probhi[2] - problo[2]);
        }

        if (perturb_vel) {
            const amrex::Real xl = x - problo[0];
            const amrex::Real yl = y - problo[1];
            const amrex::Real zl = z / ref_height;
            const amrex::Real damp = std::exp(-0.5 * zl * zl);

            velocity(i, j, k, 0) += ufac * damp * z * std::cos(aval * yl);
            velocity(i, j, k, 1) += vfac * damp * z * std::cos(bval * xl);
        }
    });
}

void ABLFieldInit::perturb_temperature(
    const int lev,
    const amrex::Geometry& geom,
    // cppcheck-suppress constParameter
    Field& temperature) const
{
    /** Perturbations for the temperature field
     *
     */

    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    auto& theta_fab = temperature(lev);
    const auto theta_cutoff_height = m_theta_cutoff_height;
    const auto theta_gauss_mean = m_theta_gauss_mean;
    const auto theta_gauss_var = m_theta_gauss_var;
    const auto deltaT = m_deltaT;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(theta_fab, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const auto& bx = mfi.tilebox();
        const auto& theta = theta_fab.array(mfi);

        amrex::ParallelForRNG(
            bx, [=] AMREX_GPU_DEVICE(
                    int i, int j, int k,
                    const amrex::RandomEngine& engine) noexcept {
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                if (z < theta_cutoff_height) {
                    theta(i, j, k) =
                        deltaT * amrex::RandomNormal(
                                     theta_gauss_mean, theta_gauss_var, engine);
                }
            });
    }
}

//! Initialize sfs tke field at the beginning of the simulation
void ABLFieldInit::init_tke(
    const amrex::Geometry& geom, amrex::MultiFab& tke_fab) const
{
    tke_fab.setVal(m_tke_init, 1);

    if (!m_tke_init_profile) {
        return;
    }

    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    const auto tke_cutoff_height = m_tke_cutoff_height;
    const auto tke_init_factor = m_tke_init_factor;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(tke_fab, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const auto& bx = mfi.growntilebox(1);
        const auto& tke = tke_fab.array(mfi);

        // Profile definition from Beare et al. (2006)
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                if (z < tke_cutoff_height) {
                    tke(i, j, k) = tke_init_factor *
                                   std::pow(1. - z / tke_cutoff_height, 3.0);
                } else {
                    tke(i, j, k) = 0.;
                }
            });
    }
}

} // namespace amr_wind
