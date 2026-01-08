#include "amr-wind/equation_systems/icns/source_terms/ABLForcing.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "amr-wind/utilities/trig_ops.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"

namespace amr_wind::pde::icns {

ABLForcing::ABLForcing(const CFDSim& sim)
    : m_time(sim.time()), m_mesh(sim.mesh())
{
    const auto& abl = sim.physics_manager().get<amr_wind::ABL>();
    abl.register_forcing_term(this);
    abl.abl_statistics().register_forcing_term(this);

    amrex::ParmParse pp_abl(identifier());
    // TODO: Allow forcing at multiple heights
    pp_abl.get("abl_forcing_height", m_forcing_height);
    amrex::ParmParse pp_incflo("incflo");

    pp_abl.query("velocity_timetable", m_vel_timetable);
    if (!m_vel_timetable.empty()) {
        std::ifstream ifh(m_vel_timetable, std::ios::in);
        if (!ifh.good()) {
            amrex::Abort(
                "Cannot find ABLForcing velocity_timetable file: " +
                m_vel_timetable);
        }
        amrex::Real data_time;
        amrex::Real data_speed;
        amrex::Real data_deg;
        ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        while (ifh >> data_time) {
            ifh >> data_speed >> data_deg;
            amrex::Real data_rad = utils::radians(data_deg);
            m_time_table.push_back(data_time);
            m_speed_table.push_back(data_speed);
            m_direction_table.push_back(data_rad);
        }
    } else {
        pp_incflo.getarr("velocity", m_target_vel);
    }

    m_write_force_timetable = pp_abl.contains("forcing_timetable_output_file");
    if (m_write_force_timetable) {
        amrex::Print() << "Here!!!!" << std::endl;
        pp_abl.get("forcing_timetable_output_file", m_force_timetable);
        pp_abl.query("forcing_timetable_frequency", m_force_outfreq);
        pp_abl.query("forcing_timetable_start_time", m_force_outstart);
        if (amrex::ParallelDescriptor::IOProcessor()) {
            amrex::Print() << m_force_timetable << std::endl;
            std::ofstream outfile;
            outfile.open(m_force_timetable, std::ios::out);
            outfile << "time\tfx\tfy\tfz\tUgx\tUgy\tzi\n";
        }
    }

    // If free-atmosphere damping is used, read these inputs.
    pp_abl.query("free_atmosphere_damping", m_fa_damping);
    if (m_fa_damping) {
        pp_abl.query("detect_free_atmosphere_height", m_fa_detect_height);
        if (!m_fa_detect_height) {
            pp_abl.get("free_atmosphere_height", m_fa_height);
        }
        pp_abl.query("free_atmosphere_damping_start_time", m_fa_time_start);
        pp_abl.query("free_atmosphere_damping_end_time", m_fa_time_end);
        pp_abl.query("free_atmosphere_damping_time_scale", m_fa_tau);

        amrex::ParmParse pp_coriolis("CoriolisForcing");
        amrex::Real rot_time_period = 86164.091;
        pp_coriolis.query("rotational_time_period", rot_time_period);

        amrex::Real omega = utils::two_pi() / rot_time_period;

        amrex::Real latitude = 90.0;
        pp_coriolis.query("latitude", latitude);
        latitude = utils::radians(latitude);
        amrex::Real sinphi = std::sin(latitude);

        m_coriolis_factor = 2.0 * omega * sinphi;
        amrex::Print() << "m_coriolis_factor = " << m_coriolis_factor << " " << omega << " " << sinphi << std::endl;
    }

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        m_mean_vel[i] = m_target_vel[i];
    }

    // Set up relaxation toward 0 forcing near the air-water interface
    if (sim.repo().field_exists("vof")) {
        // If vof exists, get multiphase physics
        const auto& mphase = sim.physics_manager().get<amr_wind::MultiPhase>();
        // Retrieve interface position
        m_water_level = mphase.water_level();
        // Confirm that water level will be used
        m_use_phase_ramp = true;
        // Parse for thickness of ramping function
        pp_abl.get("abl_forcing_off_height", m_forcing_mphase0);
        pp_abl.get("abl_forcing_ramp_height", m_forcing_mphase1);
        // Store reference to vof field
        m_vof = &sim.repo().get_field("vof");
        // Parse for number of cells in band
        pp_abl.query("abl_forcing_band", m_n_band);
    } else {
        // Point to something, will not be used
        m_vof = &sim.repo().get_field("velocity");
    }

    mean_velocity_init(abl.abl_statistics().vel_profile_coarse());
}

ABLForcing::~ABLForcing() = default;

void ABLForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{

    const amrex::Real dudt = m_abl_forcing[0];
    const amrex::Real dvdt = m_abl_forcing[1];

    const bool fa_damping = m_fa_damping;
    const bool fa_detect_height = m_fa_detect_height;
    const amrex::Real fa_height = fa_detect_height ? m_zi : m_fa_height;
    const amrex::Real fa_time_start = m_fa_time_start;
    const amrex::Real fa_time_end = m_fa_time_end;
    const amrex::Real fa_tau = m_fa_tau;
    const amrex::Real fa_u = dvdt / m_coriolis_factor;
    const amrex::Real fa_v = -dudt / m_coriolis_factor;

    const auto& current_time = m_time.current_time();
    const auto& new_time = m_time.new_time();
    const auto& nph_time = 0.5 * (current_time + new_time);

    const bool ph_ramp = m_use_phase_ramp;
    const int n_band = m_n_band;
    const amrex::Real wlev = m_water_level;
    const amrex::Real wrht0 = m_forcing_mphase0;
    const amrex::Real wrht1 = m_forcing_mphase1;
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();

    const int idir = m_axis;
    const amrex::Real* heights = m_vel_ht.data();
    const amrex::Real* heights_end = m_vel_ht.end();
    const amrex::Real* vals = m_vel_vals.data();

    const auto& vof = (*m_vof)(lev).const_array(mfi);

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::IntVect iv(i, j, k);
        const amrex::Real z = problo[idir] + (iv[idir] + 0.5) * dx[idir];

        const amrex::Real umean =
            amr_wind::interp::linear(heights, heights_end, vals, z, 3, 0);
        const amrex::Real vmean =
            amr_wind::interp::linear(heights, heights_end, vals, z, 3, 1);

        // This part deals with air and water phase and only masks ABL forcing
        // on the water phase.
        amrex::Real fac = 1.0;
        if (ph_ramp) {
            if (z - wlev < wrht0 + wrht1) {
                if (z - wlev < wrht0) {
                    // Apply no forcing within first interval
                    fac = 0.0;
                } else {
                    // Ramp from 0 to 1 over second interval
                    fac =
                        0.5 - 0.5 * std::cos(M_PI * (z - wlev - wrht0) / wrht1);
                }
            }
            // Check for presence of liquid (like a droplet)
            // - interface_band checks for closeness to interface
            // - need to also check for within liquid
            if (multiphase::interface_band(i, j, k, vof, n_band) ||
                vof(i, j, k) > 1.0 - 1e-12) {
                // Turn off forcing
                fac = 0.0;
            }
        }

        // This part applies the free atmosphere damping, which nudges
        // the solution in the free atmosphere toward geostrophic as
        // computed by the wind speed controller. (i.e., the wind speed
        // controller computes a forcing term that is really a driving
        // pressure gradient, so there is a corresponding geostrophic
        // wind. Nudge the free atmosphere toward that geostrophic wind.)
        amrex::Real src_fa_damping_x = 0.0;
        amrex::Real src_fa_damping_y = 0.0;
        if (fa_damping &&
            ((nph_time >= fa_time_start) && (nph_time <= fa_time_end)) &&
            (z >= fa_height)) {
            src_fa_damping_x = (1.0 / fa_tau) * (fa_u - umean);
            src_fa_damping_y = (1.0 / fa_tau) * (fa_v - vmean);
        }

        // Sum up the source term
        src_term(i, j, k, 0) += (fac * dudt) + src_fa_damping_x;
        src_term(i, j, k, 1) += (fac * dvdt) + src_fa_damping_y;
        src_term(i, j, k, 2) += 0.0;
    });
}

void ABLForcing::mean_velocity_init(const VelPlaneAveraging& vavg)
{
    m_axis = vavg.axis();

    // The implementation depends the assumption that the ABL statistics class
    // computes statistics at the cell-centeres only on level 0. If this
    // assumption changes in future, the implementation will break... so put in
    // a check here to catch this.
    AMREX_ALWAYS_ASSERT(
        m_mesh.Geom(0).Domain().length(m_axis) ==
        static_cast<int>(vavg.line_centroids().size()));

    m_vel_ht.resize(vavg.line_centroids().size());
    m_vel_vals.resize(vavg.line_average().size());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, vavg.line_centroids().begin(),
        vavg.line_centroids().end(), m_vel_ht.begin());

    mean_velocity_update(vavg);
}

void ABLForcing::mean_velocity_update(const VelPlaneAveraging& vavg)
{
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, vavg.line_average().begin(),
        vavg.line_average().end(), m_vel_vals.begin());
}

} // namespace amr_wind::pde::icns
