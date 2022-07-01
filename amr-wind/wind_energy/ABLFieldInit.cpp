#include <cmath>

#include "amr-wind/wind_energy/ABLFieldInit.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

ABLFieldInit::ABLFieldInit(CFDSim& sim)
    : m_repo(sim.repo()), m_mesh_mapping(sim.has_mesh_mapping())
{
    amrex::ParmParse pp_abl("ABL");
    amrex::ParmParse pp_ablf("ABLForcing");

    // Temperature variation as a function of height
    pp_abl.getarr("temperature_heights", m_theta_heights);
    pp_abl.getarr("temperature_values", m_theta_values);

    AMREX_ALWAYS_ASSERT(m_theta_heights.size() == m_theta_values.size());
    int num_theta_values = m_theta_heights.size();

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

    // Initial profile type: const, simple
    // TODO: Add M-O and geostrophic forcing-consistent profiles 
    pp_abl.query("init_type", m_init_type);
    pp_abl.query("kappa", m_kappa);
    pp_abl.query("surface_roughness_z0", m_rough_z0);
    pp_ablf.query("abl_forcing_height", m_force_height);
    pp_abl.query("init_tke_const", m_tke_init);
    pp_abl.query("init_sdr_const", m_sdr_init);
    pp_abl.query("init_eps_const", m_eps_init);
    pp_abl.query("init_profile_ustar", m_init_prof_ustar);
    pp_abl.query("init_profile_c1", m_init_prof_c1);
    pp_abl.query("init_profile_c2", m_init_prof_c2);

    // Get beta-star for use in initial conditions
    amrex::ParmParse pp_turb("turbulence");
    pp_turb.query("model", m_turb_model_name);
    const std::string turb_coeffs_dict = m_turb_model_name + "_coeffs";
    amrex::ParmParse pp_turb_coeffs(turb_coeffs_dict);
    pp_turb.query("beta_star", m_beta_star);

    // TODO: Modify this to accept velocity as a function of height
    // Extract velocity field from incflo
    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.get("density", m_rho);
    pp_incflo.getarr("velocity", m_vel);

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
    const int level,
    const amrex::MFIter& mfi,
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
    const amrex::Real rho_init = m_rho;
    const amrex::Real umean = m_vel[0];
    const amrex::Real vmean = m_vel[1];
    const amrex::Real wmean = m_vel[2];
    const amrex::Real aval = m_Uperiods * 2.0 * pi / (probhi[1] - problo[1]);
    const amrex::Real bval = m_Vperiods * 2.0 * pi / (probhi[0] - problo[0]);
    const amrex::Real ufac = m_deltaU * std::exp(0.5) / m_ref_height;
    const amrex::Real vfac = m_deltaV * std::exp(0.5) / m_ref_height;
    const amrex::Real ref_height = m_ref_height;

    const int ntvals = m_theta_heights.size();
    const amrex::Real* th = m_thht_d.data();
    const amrex::Real* tv = m_thvv_d.data();

    const amrex::Real kappa = m_kappa;
    const amrex::Real z_fh = m_force_height;
    const amrex::Real z_0 = m_rough_z0;

    // Longitudinal velocity angle (positive from X-dir)
    const amrex::Real vel_angle = std::atan(vmean/umean);

    // Longitudinal Velocity Magnitude
    const amrex::Real vel_long = std::sqrt(std::pow(umean,2.0) + std::pow(vmean, 2.0));

    // Calc ustar
    const amrex::Real ustar = vel_long*kappa/std::log((z_fh + z_0)/z_0);
    amrex::Print() << "Calculated u*: " << ustar << std::endl;

    // Handle mesh mapping coordinates
    Field const* nu_coord_cc =
        m_mesh_mapping ? &(m_repo.get_field("non_uniform_coord_cc")) : nullptr;
    amrex::Array4<amrex::Real const> nu_cc =
        m_mesh_mapping ? ((*nu_coord_cc)(level).array(mfi))
                       : amrex::Array4<amrex::Real const>();

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
	    const amrex::Real x = 
            m_mesh_mapping ? nu_cc(i, j, k, 0) : problo[0] + (i + 0.5) * dx[0];
	    const amrex::Real y = 
            m_mesh_mapping ? nu_cc(i, j, k, 1) : problo[1] + (j + 0.5) * dx[1];
	    const amrex::Real z = 
            m_mesh_mapping ? nu_cc(i, j, k, 2) : problo[2] + (k + 0.5) * dx[2];

        density(i, j, k) = rho_init;

        // Mean velocity field
        if(m_init_type == "const"){
            velocity(i, j, k, 0) = umean;
            velocity(i, j, k, 1) = vmean;
            velocity(i, j, k, 2) = wmean;
        }else{
            velocity(i, j, k, 0) = ustar/kappa*std::log((z+z_0)/z_0)*std::cos(vel_angle);
            velocity(i, j, k, 1) = ustar/kappa*std::log((z+z_0)/z_0)*std::sin(vel_angle);
            velocity(i, j, k, 2) = wmean;
        }

        amrex::Real theta = tv[0];
        for (int iz = 0; iz < ntvals - 1; ++iz) {
            if ((z > th[iz]) && (z <= th[iz + 1])) {
                const amrex::Real slope =
                    (tv[iz + 1] - tv[iz]) / (th[iz + 1] - th[iz]);
                theta = tv[iz] + (z - th[iz]) * slope;
            }
        }

        temperature(i, j, k, 0) += theta;

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
    /** Perturbations for the temperature field is adapted from the following
     * paper:
     *
     *  D. Munoz-Esparza, B. Kosovic, J. van Beeck, J. D. Mirocha, A stocastic
     *  perturbation method to generate inflow turbulence in large-eddy
     *  simulation models: Application to neutrally stratified atmospheric
     *  boundary layers. Physics of Fluids, Vol. 27, 2015.
     *
     */

    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    auto& theta_fab = temperature(lev);
    const auto theta_cutoff_height = m_theta_cutoff_height;
    const auto theta_gauss_mean = m_theta_gauss_mean;
    const auto theta_gauss_var = m_theta_gauss_var;
    const auto deltaT = m_deltaT;
    Field const* nu_coord_cc =
        m_mesh_mapping ? &(m_repo.get_field("non_uniform_coord_cc")) : nullptr;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(theta_fab, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const auto& bx = mfi.tilebox();
        const auto& theta = theta_fab.array(mfi);
        amrex::Array4<amrex::Real const> nu_cc =
            m_mesh_mapping ? ((*nu_coord_cc)(lev).array(mfi))
                           : amrex::Array4<amrex::Real const>();
        amrex::ParallelForRNG(
            bx, [=] AMREX_GPU_DEVICE(
                    int i, int j, int k,
                    const amrex::RandomEngine& engine) noexcept {
                const amrex::Real z = m_mesh_mapping 
                                          ? nu_cc(i, j, k, 2)
                                          : problo[2] + (k + 0.5) * dx[2];

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
    const int level,
    const amrex::MFIter& mfi,
    const amrex::Box& vbx,
    const amrex::Geometry& geom,
    const amrex::Array4<amrex::Real>& tke) const
{
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    const auto& probhi = geom.ProbHiArray();

    // Handle mesh mapping coordinates
    Field const* nu_coord_cc =
        m_mesh_mapping ? &(m_repo.get_field("non_uniform_coord_cc")) : nullptr;
    amrex::Array4<amrex::Real const> nu_cc =
        m_mesh_mapping ? ((*nu_coord_cc)(level).array(mfi))
                    : amrex::Array4<amrex::Real const>();

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real x = 
            m_mesh_mapping ? nu_cc(i, j, k, 0) : problo[0] + (i + 0.5) * dx[0];
        const amrex::Real y = 
            m_mesh_mapping ? nu_cc(i, j, k, 1) : problo[1] + (j + 0.5) * dx[1];
        const amrex::Real z = 
            m_mesh_mapping ? nu_cc(i, j, k, 2) : problo[2] + (k + 0.5) * dx[2];

        if(m_init_type == "const"){
            tke(i, j, k, 0) = m_tke_init;
        }else{
            tke(i, j, k, 0) = m_init_prof_c1*std::log(z + m_rough_z0) + m_init_prof_c2;
        }

    });
}

//! Initialize SDR field at the beginning of the simulation
//! (applicable to K-Omega SST model)
void ABLFieldInit::init_sdr(
    const int level,
    const amrex::MFIter& mfi,
    const amrex::Box& vbx,
    const amrex::Geometry& geom,
    const amrex::Array4<amrex::Real>& sdr) const
{
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    const auto& probhi = geom.ProbHiArray();

    // Longitudinal velocity angle (positive from X-dir)
    const amrex::Real vel_angle = std::atan(m_vel[1]/m_vel[0]);

    // Longitudinal Velocity Magnitude
    const amrex::Real vel_long = std::sqrt(std::pow(m_vel[0],2.0) + std::pow(m_vel[1], 2.0));

    // Calc ustar
    const amrex::Real ustar = vel_long*m_kappa/std::log((m_force_height + m_rough_z0)/m_rough_z0);

    // Handle mesh mapping coordinates
    Field const* nu_coord_cc =
        m_mesh_mapping ? &(m_repo.get_field("non_uniform_coord_cc")) : nullptr;
    amrex::Array4<amrex::Real const> nu_cc =
        m_mesh_mapping ? ((*nu_coord_cc)(level).array(mfi))
                    : amrex::Array4<amrex::Real const>();

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real x = 
            m_mesh_mapping ? nu_cc(i, j, k, 0) : problo[0] + (i + 0.5) * dx[0];
        const amrex::Real y = 
            m_mesh_mapping ? nu_cc(i, j, k, 1) : problo[1] + (j + 0.5) * dx[1];
        const amrex::Real z = 
            m_mesh_mapping ? nu_cc(i, j, k, 2) : problo[2] + (k + 0.5) * dx[2];

        if(m_init_type == "const"){
            sdr(i, j, k, 0) = m_sdr_init;
        }else{
            const amrex::Real initk = m_init_prof_c1*std::log(z + m_rough_z0) + m_init_prof_c2;
            sdr(i, j, k, 0) = std::pow(ustar, 3.0)/(m_beta_star*initk*m_kappa*(z + m_rough_z0));
        }

    });
}

//! Initialize epsilon field at the beginning of the simulation
//! (applicable to K-Epsilon model)
void ABLFieldInit::init_eps(
    const int level,
    const amrex::MFIter& mfi,
    const amrex::Box& vbx,
    const amrex::Geometry& geom,
    const amrex::Array4<amrex::Real>& eps) const
{
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    const auto& probhi = geom.ProbHiArray();

    // Longitudinal velocity angle (positive from X-dir)
    const amrex::Real vel_angle = std::atan(m_vel[1]/m_vel[0]);

    // Longitudinal Velocity Magnitude
    const amrex::Real vel_long = std::sqrt(std::pow(m_vel[0],2.0) + std::pow(m_vel[1], 2.0));

    // Calc ustar
    const amrex::Real ustar = vel_long*m_kappa/std::log((m_force_height + m_rough_z0)/m_rough_z0);

    // Handle mesh mapping coordinates
    Field const* nu_coord_cc =
        m_mesh_mapping ? &(m_repo.get_field("non_uniform_coord_cc")) : nullptr;
    amrex::Array4<amrex::Real const> nu_cc =
        m_mesh_mapping ? ((*nu_coord_cc)(level).array(mfi))
                    : amrex::Array4<amrex::Real const>();

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real x = 
            m_mesh_mapping ? nu_cc(i, j, k, 0) : problo[0] + (i + 0.5) * dx[0];
        const amrex::Real y = 
            m_mesh_mapping ? nu_cc(i, j, k, 1) : problo[1] + (j + 0.5) * dx[1];
        const amrex::Real z = 
            m_mesh_mapping ? nu_cc(i, j, k, 2) : problo[2] + (k + 0.5) * dx[2];

        if(m_init_type == "const"){
            eps(i, j, k, 0) = m_sdr_init;
        }else{
            eps(i, j, k, 0) = std::pow(ustar, 3.0)/(m_kappa*(z + m_rough_z0));
        }

    });
}


} // namespace amr_wind
