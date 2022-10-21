#include <cmath>

#include "amr-wind/wind_energy/ABLFieldInit.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"

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

    // Use input from netcdf file
    pp_abl.query("initial_condition_input_file", m_ic_input);

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
    const amrex::Real rho_init = m_rho;

    const amrex::Real umean =
        !m_vel_timetable.empty() ? m_vel_speed * std::cos(m_vel_dir) : m_vel[0];
    const amrex::Real vmean =
        !m_vel_timetable.empty() ? m_vel_speed * std::sin(m_vel_dir) : m_vel[1];
    const amrex::Real wmean = !m_vel_timetable.empty() ? 0.0 : m_vel[2];

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

        if (perturb_vel) {
            const amrex::Real xl = x - problo[0];
            const amrex::Real yl = y - problo[1];
            const amrex::Real zl = z / ref_height;
            const amrex::Real damp = std::exp(-0.5 * zl * zl);

            velocity(i, j, k, 0) += ufac * damp * z * std::cos(aval * yl);
            velocity(i, j, k, 1) += vfac * damp * z * std::cos(bval * xl);
        }
    });

#ifdef AMR_WIND_USE_NETCDF
    // Load the netcdf file with data if specified in the inputs
    if (!m_ic_input.empty()) {

        // Open the netcdf input file
        // This file should have the same dimensions as the simulation
        auto ncf = ncutils::NCFile::open(m_ic_input, NC_NOWRITE);

        // Ensure that the input dimensions match the coarsest grid size
        const auto& domain = geom.Domain();

        // The indices that determine the start and end points of the i, j, k
        // arrays. The max and min are there to ensure that the points form
        // ghost cells are not used
        auto i0 = std::max(vbx.smallEnd(0), domain.smallEnd(0));
        auto i1 = std::min(vbx.bigEnd(0), domain.bigEnd(0));

        auto j0 = std::max(vbx.smallEnd(1), domain.smallEnd(1));
        auto j1 = std::min(vbx.bigEnd(1), domain.bigEnd(1));

        auto k0 = std::max(vbx.smallEnd(2), domain.smallEnd(2));
        auto k1 = std::min(vbx.bigEnd(2), domain.bigEnd(2));
        // std::cout << "k vals" <<  k0 << " " << k1 << std::endl;
        // The x, y and z velocity components (u, v, w)
        auto uvel = ncf.var("uvel");
        auto vvel = ncf.var("vvel");

        // Loop through all points in the domain and set velocities to values
        // from the input file
        // start is the first index from where to read data
        std::vector<size_t> start{
            {static_cast<size_t>(i0), static_cast<size_t>(j0),
             static_cast<size_t>(k0)}};
        // count is the total number of elements to read in each direction
        std::vector<size_t> count{
            {static_cast<size_t>(i1 - i0 + 1), static_cast<size_t>(j1 - j0 + 1),
             static_cast<size_t>(k1 - k0 + 1)}};

        // Vector to store the 3d data into a single array
        amrex::Vector<double> uvel2;
        amrex::Vector<double> vvel2;

        // Set the size of the arrays to the total number of points in this
        // processor
        uvel2.resize(count[0] * count[1] * count[2]);
        vvel2.resize(count[0] * count[1] * count[2]);

        // Read the velocity components u and v
        uvel.get(uvel2.data(), start, count);
        vvel.get(vvel2.data(), start, count);

        // Amrex parallel for to assign the velocity at each point
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // The counter to go from 3d to 1d vector
                // auto idx = i + j * count[1] + k * count[2] * count[1];
                // auto idx = (i - i0) + (j - j0) * count[1] + (k - k0) *
                // count[2] * count[1];
                // auto idx = i * count[0] * count[1] +  j * count[1]  + k;
                auto idx = (i - i0) * count[2] * count[1] +
                           (j - j0) * count[2] + (k - k0);
                // auto idx = (i-i0) +  (j-j0) * count[0]  + (k-k0) * count[0] *
                // count[1];
                velocity(i, j, k, 0) = uvel2.data()[idx];
                velocity(i, j, k, 1) = vvel2.data()[idx];
            });
        // Close the netcdf file
        ncf.close();
    }
// Make sure to call fill patch *******
#endif
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
    const amrex::Geometry& /* geom */, amrex::MultiFab& tke) const
{
    tke.setVal(m_tke_init, 1);
}

} // namespace amr_wind
