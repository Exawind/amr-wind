#include "amr-wind/utilities/sampling/RadarSampler.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/utilities/trig_ops.H"

#include "AMReX_ParmParse.H"

#define PI 3.14159265

namespace amr_wind::sampling {

RadarSampler::RadarSampler(const CFDSim& sim) : m_sim(sim) {}

RadarSampler::~RadarSampler() = default;

void RadarSampler::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);

    // This is the origin of the scan (x, y, z) [m]
    pp.getarr("origin", m_start);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_start.size()) == AMREX_SPACEDIM);

    // Sampling Frequency
    pp.get("sampling_frequency", m_sample_freq);

    // Sampling Frequency
    pp.get("device_sampling_frequency", m_radar_sample_freq);

    // The number of points in single beam
    pp.get("num_points", m_npts);

    // Half angle of the cone triangle [deg]
    pp.get("radar_cone_angle", m_cone_angle);

    // Quadrature type ("truncated_normal_halfpower" supported)
    pp.get("radar_quadrature_type", m_radar_quad_type);

    // Number of radials in the cone [-]
    pp.get("radar_npts_azimuth", m_npts_azimuth);

    // The length of the beam [m]
    pp.get("radar_beam_length", m_beam_length);

    // Angular speed of the beam [deg/s]
    pp.get("angular_speed", m_angular_speed);

    // Unit vector, horizontal center-axis of the beam [-]
    pp.getarr("axis", m_axis);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_axis.size()) == AMREX_SPACEDIM);

    // Unit vertical axis [-]
    pp.getarr("vertical_unit_dir", m_vertical);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_vertical.size()) == AMREX_SPACEDIM);

    // Sweep angle of the beam [deg]
    pp.get("sweep_angle", m_sweep_angle);

    // Reset time of the beam at elevation change [s]
    pp.get("reset_time", m_reset_time);

    // Elevation angles [deg]
    pp.getarr("elevation_angles", m_elevation_angles);

    // Print radar debug info
    pp.query("debug_print", m_debug_print);

    // Output cone points
    pp.query("output_cone_points", m_output_cone_points);

    // Set initial m_end but isn't used for anything
    m_end.resize(AMREX_SPACEDIM);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        m_end[d] = m_axis[d] * m_beam_length + m_start[d];
    }

    // Set initial periodic time to half-way point of forward sweep
    // But this gets overwritten in update_sampling_locations
    m_periodic_time = 0.5 * m_sweep_angle / m_angular_speed;
    m_current_phase = determine_operation_phase();

    m_cone_size = num_points_cone();

    // Calculate total number of necessary beams per timstep
    // Due to subsampling freq
    amrex::Real dt_sim = m_sim.time().delta_t();
    amrex::Real dt_sample = 1.0 / m_sample_freq;
    double timestep_sample_ratio = dt_sim / dt_sample;
    m_ntotal = int(std::ceil(timestep_sample_ratio));

    AMREX_ALWAYS_ASSERT(timestep_sample_ratio >= 1.0);

    m_prior_cones.resize(num_points_scan());
    m_current_cones.resize(num_points_scan());

    RadarSampler::new_cone();

    check_bounds();
}

void RadarSampler::check_bounds()
{
    const int lev = 0;
    const auto* prob_lo = m_sim.mesh().Geom(lev).ProbLo();
    const auto* prob_hi = m_sim.mesh().Geom(lev).ProbHi();

    bool all_ok = true;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (m_start[d] < prob_lo[d]) {
            all_ok = false;
            m_start[d] = prob_lo[d];
        }
        if (m_start[d] > prob_hi[d]) {
            all_ok = false;
            m_start[d] = prob_lo[d];
        }
        if (m_end[d] < prob_lo[d]) {
            all_ok = false;
            m_end[d] = prob_lo[d];
        }
        if (m_end[d] > prob_hi[d]) {
            all_ok = false;
            m_end[d] = prob_lo[d];
        }
    }
    if (!all_ok) {
        amrex::Print() << "WARNING: RadarSampler: Out of domain line was "
                          "truncated to match domain"
                       << std::endl;
    }
}

double RadarSampler::total_sweep_time() const
{
    return 2 * (m_sweep_angle / m_angular_speed + m_reset_time);
}

double RadarSampler::periodic_time()
{
    m_periodic_time =
        m_radar_time -
        std::floor(m_radar_time / total_sweep_time()) * total_sweep_time();
    return m_periodic_time;
}

// Single direction sweep count
int RadarSampler::sweep_count() const
{
    const double sweep_time = m_sweep_angle / m_angular_speed + m_reset_time;
    return static_cast<int>(std::floor(m_radar_time / sweep_time));
}

RadarSampler::phase RadarSampler::determine_operation_phase() const
{
    const double phase_time = m_sweep_angle / m_angular_speed;

    RadarSampler::phase current_phase;

    if (m_periodic_time < phase_time) {
        current_phase = phase::FORWARD;
    } else if (m_periodic_time < phase_time + m_reset_time) {
        current_phase = phase::FORWARD_PAUSE;
    } else if (m_periodic_time < 2 * phase_time + m_reset_time) {
        current_phase = phase::REVERSE;
    } else {
        current_phase = phase::REVERSE_PAUSE;
    }

    return current_phase;
}

double RadarSampler::determine_current_sweep_angle() const
{
    switch (determine_operation_phase()) {
    case phase::FORWARD: {
        return m_angular_speed * m_periodic_time - m_sweep_angle / 2;
    }
    case phase::FORWARD_PAUSE: {
        return m_sweep_angle / 2;
    }
    case phase::REVERSE: {
        return 3 * m_sweep_angle / 2 -
               m_angular_speed * (m_periodic_time - m_reset_time);
    }
    default: {
        return -m_sweep_angle / 2;
    }
    }
}

bool RadarSampler::update_sampling_locations()
{
    amrex::Real time = m_sim.time().current_time();
    amrex::Real start_time = m_sim.time().start_time();
    amrex::Real dt_sim = m_sim.time().delta_t();
    amrex::Real dt_sample = 1.0 / m_sample_freq;
    amrex::Real st_diff = time - start_time;

    if (st_diff < 1.0e-10) {
        m_radar_time = time;
    }

    amrex::Real ts_diff = time - m_radar_time;

    // Correction for time mismatch
    int time_corr = (ts_diff > dt_sample) ? int(ts_diff / dt_sample) : 0;

    m_ns = int(dt_sim / dt_sample) + time_corr;

    m_los_unit.resize(m_ntotal * num_points_quad());
    m_los_proj.resize(m_ntotal * num_points_quad());

    if (m_debug_print) {
        amrex::Print() << "-------------------------" << std::endl
                       << "ts_diff: " << ts_diff << "\tm_ns: " << m_ns
                       << "\tm_ntotal: " << m_ntotal
                       << "\tRadar Time: " << m_radar_time
                       << "\tAMR Time: " << time << std::endl
                       << "-------------------------" << std::endl;
    }

    long conetipbegin = m_cone_size - 1 - num_points_quad();
    long conetipend = m_cone_size;

    // Loop for oversampling
    for (int k = 0; k < m_ntotal; k++) {
        if (k < m_ns) {
            double per_time = periodic_time();
            double sweep_angle = determine_current_sweep_angle();
            double elevation_angle = m_elevation_angles.at(
                sweep_count() % m_elevation_angles.size());

            if (m_debug_print) {
                amrex::Print() << "-------------------------" << std::endl
                               << "Total Sweep: " << total_sweep_time() << "\t"
                               << "Periodic Time: " << per_time << "\t"
                               << "Radar Time: " << m_radar_time << "\t"
                               << "S Angle: " << sweep_angle << "\t"
                               << "E Angle: " << elevation_angle << "\t"
                               << "Sim Time: " << time << std::endl
                               << "-------------------------" << std::endl;
            }

            vs::Vector vertical_ref(
                m_vertical[0], m_vertical[1], m_vertical[2]);
            vs::Vector radar_ref(m_axis[0], m_axis[1], m_axis[2]);

            // Assume vertical_ref and radar_ref are normal to each other and
            // unit
            vs::Vector swept_axis(sampling_utils::rotate_euler_vector(
                vertical_ref, sweep_angle, radar_ref));
            vs::Vector elevation_axis(-vertical_ref ^ swept_axis);

            long cq_idx = k * num_points_quad();

            // Add rotated cone to current cones
            for (int i = 0; i < m_cone_size; ++i) {
                vs::Vector temp_point(
                    m_initial_cone[i][0] - m_start[0],
                    m_initial_cone[i][1] - m_start[1],
                    m_initial_cone[i][2] - m_start[2]);
                vs::Vector swept_point(sampling_utils::rotate_euler_vector(
                    vertical_ref, sweep_angle, temp_point));
                vs::Vector rotated_point(sampling_utils::rotate_euler_vector(
                    elevation_axis, elevation_angle, swept_point));
                long point_index = i + k * m_cone_size;
                m_current_cones[point_index][0] = rotated_point[0] + m_start[0];
                m_current_cones[point_index][1] = rotated_point[1] + m_start[1];
                m_current_cones[point_index][2] = rotated_point[2] + m_start[2];

                // Add a single cap to help calc line of sight
                if (i > conetipbegin && i < conetipend) {
                    vs::Vector unit_cone_point(rotated_point.normalize());
                    vs::Tensor unit_proj_mat =
                        sampling_utils::unit_projection_matrix(unit_cone_point);
                    m_los_unit[cq_idx] = unit_cone_point;
                    m_los_proj[cq_idx] = unit_proj_mat;
                    cq_idx++;
                }
            }

            m_radar_time += dt_sample;

        } else {
            // This cone falls outside time bounds
            // For this timestep
            long cq_idx = k * num_points_quad();

            for (int i = 0; i < m_cone_size; ++i) {
                long point_index = i + k * m_cone_size;
                m_current_cones[point_index][0] = m_fill_val;
                m_current_cones[point_index][1] = m_fill_val;
                m_current_cones[point_index][2] = m_fill_val;

                // Create zero projection matrix
                if (i > conetipbegin && i < conetipend) {
                    vs::Vector unit_zero_point(0.0, 0.0, 0.0);
                    vs::Tensor zero_mat =
                        sampling_utils::unit_projection_matrix(unit_zero_point);
                    m_los_unit[cq_idx] = unit_zero_point;
                    m_los_proj[cq_idx] = zero_mat;
                    cq_idx++;
                }
            }
        }
    }

    m_radar_iter++;

    return true;
}

void RadarSampler::new_cone()
{
    vs::Vector canon_vector(0.0, 0.0, 1.0);
    vs::Vector origin_vector(0.0, 0.0, 0.0);
    vs::Vector init_axis(m_axis[0], m_axis[1], m_axis[2]);

    m_initial_cone.resize(num_points_cone());
    m_rays.resize(num_points_quad());
    m_weights.resize(num_points_quad());

    sampling_utils::spherical_cap_truncated_normal(
        utils::radians(m_cone_angle), m_ntheta,
        sampling_utils::NormalRule::HALFPOWER, m_rays, m_weights);

    long nquad = num_points_quad();

    if (m_debug_print) {
        for (int j = 0; j < nquad; ++j) {
            amrex::Print() << "Rays x: " << m_rays[j][0] << "\t"
                           << "Rays y: " << m_rays[j][1] << "\t"
                           << "Rays z: " << m_rays[j][2] << "\t" << std::endl;
        }
    }

    amrex::Array<amrex::Real, AMREX_SPACEDIM> dx;

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        dx[d] = (canon_vector[d] - origin_vector[d]) / m_npts;
    }

    if (m_debug_print) {
        amrex::Print() << "-------------------------" << std::endl
                       << "Beam Length: " << m_beam_length
                       << "\tnquad: " << nquad << "\tdx[2]: " << dx[2]
                       << "\tm_npts: " << m_npts << "\tntheta: " << m_ntheta
                       << "\tnphi: " << m_nphi << std::endl
                       << "-------------------------" << std::endl;
    }

    for (int n = 0; n < m_npts; ++n) {
        for (int j = 0; j < nquad; ++j) {
            long pt_idx = j + n * nquad;
            const auto radius = dx[2] * n;
            m_initial_cone[pt_idx][0] = radius * m_rays[j][0];
            m_initial_cone[pt_idx][1] = radius * m_rays[j][1];
            m_initial_cone[pt_idx][2] = radius * m_rays[j][2];
        }
    }

    vs::Vector tc_axis = canon_vector ^ init_axis;
    amrex::Real tc_angle = std::acos(init_axis & canon_vector) * 180.0 / PI;

    // Initial cone points along input-file-given axis and is full size
    for (int i = 0; i < num_points_cone(); ++i) {
        vs::Vector temp_loc(
            m_initial_cone[i][0], m_initial_cone[i][1], m_initial_cone[i][2]);
        vs::Vector new_rot(
            sampling_utils::rotate_euler_vector(tc_axis, tc_angle, temp_loc));
        m_initial_cone[i][0] = new_rot[0] * m_beam_length + m_start[0];
        m_initial_cone[i][1] = new_rot[1] * m_beam_length + m_start[1];
        m_initial_cone[i][2] = new_rot[2] * m_beam_length + m_start[2];
    }
}

void RadarSampler::calc_lineofsight_velocity(
    const std::vector<std::vector<double>>& velocity_raw, const int interp_idx)
{
    AMREX_ALWAYS_ASSERT(
        static_cast<int>(velocity_raw.size()) == num_points_scan());

    m_los_velocity_next.resize(num_output_points());
    m_los_velocity.resize(num_output_points());
    m_los_velocity_prior.resize(num_output_points());

    std::vector<double> los_temp(num_points_scan());

    for (int k = 0; k < m_ntotal; k++) {
        for (int i = 0; i < m_cone_size; ++i) {
            long p_idx = i + k * m_cone_size;
            long cq_idx = i % num_points_quad();
            vs::Vector temp_vel(
                velocity_raw[p_idx][0], velocity_raw[p_idx][1],
                velocity_raw[p_idx][2]);
            vs::Vector los_vel_vector(temp_vel & m_los_proj[cq_idx]);
            los_temp[p_idx] = (los_vel_vector & m_los_unit[cq_idx]);
        }
    }

    for (int k = 0; k < m_ntotal; k++) {
        long a_start = k * num_points_axis();
        std::vector<double> temp_vals(
            los_temp.begin() + k * num_points_cone(),
            los_temp.begin() + (k + 1) * num_points_cone());
        if (interp_idx > 0) {
            line_average(m_weights, temp_vals, m_los_velocity_next, a_start);
        } else {
            line_average(m_weights, temp_vals, m_los_velocity, a_start);
        }
    }
}

// TODO: Fix modify_sample_data...single var output, not multiple
std::vector<double> RadarSampler::modify_sample_data(
    const std::vector<double>& sample_data, const std::string& /*unused*/)
{
    // sample_data enters this method for each sampled variable
    // there are m_ntotal steps (cones) based on sampling rate
    AMREX_ALWAYS_ASSERT(static_cast<int>(sample_data.size()) == num_points());

    const long n_cones = m_ntotal;
    std::vector<double> mod_data(num_output_points());

    for (int ic = 0; ic < n_cones; ++ic) {
        long c_start = ic * num_points_cone();
        long c_end = (ic + 1) * num_points_cone();
        long a_start = ic * num_points_axis();

        // Send a single cone to be line averaged
        std::vector<double> temp_vals(
            sample_data.begin() + c_start, sample_data.begin() + c_end);
        line_average(m_weights, temp_vals, mod_data, a_start);
    }

    return mod_data;
}

void RadarSampler::line_average(
    const std::vector<double>& weights,
    const std::vector<double>& values,
    std::vector<double>& reduced,
    long offset)
{
    const int nline = static_cast<int>(values.size() / weights.size());
    const int nquad = static_cast<int>(weights.size());
    for (int n = 0; n < nline; ++n) {
        double weight_sum = 0;
        double average = 0;
        for (int j = 0; j < nquad; ++j) {
            const int point_idx = nquad * n + j;
            weight_sum += weights[j];
            average += weights[j] * values[point_idx];
        }
        if (weight_sum > 0) {
            reduced[n + offset] = average / weight_sum;
        } else {
            reduced[n + offset] = 0;
        }
    }
}

void RadarSampler::sampling_locations(SampleLocType& locs) const
{
    locs.resize(num_points());

    for (int i = 0; i < num_points_scan(); ++i) {
        locs[i][0] = (m_radar_iter > 0) ? m_prior_cones[i][0] : m_fill_val;
        locs[i][1] = (m_radar_iter > 0) ? m_prior_cones[i][1] : m_fill_val;
        locs[i][2] = (m_radar_iter > 0) ? m_prior_cones[i][2] : m_fill_val;
    }

    for (int i = 0; i < num_points_scan(); ++i) {
        long ioff = i + num_points_scan();
        locs[ioff][0] = m_current_cones[i][0];
        locs[ioff][1] = m_current_cones[i][1];
        locs[ioff][2] = m_current_cones[i][2];
    }
}

// TODO: Create test for this calculation
void RadarSampler::cone_axis_locations(SampleLocType& axis_locs) const
{
    axis_locs.resize(num_output_points());

    for (int k = 0; k < m_ntotal; k++) {
        for (int i = 0; i < m_npts; ++i) {
            int pi = i + k * m_npts;
            long ci = i * num_points_quad() + k * num_points_cone();
            axis_locs[pi][0] = m_prior_cones[ci][0];
            axis_locs[pi][1] = m_prior_cones[ci][1];
            axis_locs[pi][2] = m_prior_cones[ci][2];
        }
    }
}

void RadarSampler::post_sample_actions()
{

    m_prior_cones.resize(num_points_scan());
    m_prior_cones = m_current_cones;

    m_los_velocity_prior.resize(num_output_points());
    m_los_velocity_prior = m_los_velocity_next;

    m_ns_prior = m_ns;
    m_ntotal_prior = m_ntotal;
    m_num_points_prior = num_points();
    m_num_output_points_prior = num_output_points();
}

#ifdef AMR_WIND_USE_NETCDF
void RadarSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
{
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("start", m_start);
    grp.put_attr("end", m_end);

    grp.def_dim("ndim", AMREX_SPACEDIM);
    grp.def_dim("num_points_cone", num_points_scan());
    grp.def_var("points", NC_DOUBLE, {"num_time_steps", "num_points", "ndim"});

    if (m_output_cone_points) {
        grp.def_var(
            "conepoints", NC_DOUBLE,
            {"num_time_steps", "num_points_cone", "ndim"});
    }
}

void RadarSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}

void RadarSampler::output_netcdf_data(
    const ncutils::NCGroup& grp, const size_t nt) const
{
    std::vector<size_t> start{nt, 0, 0};
    std::vector<size_t> count{1, 0, AMREX_SPACEDIM};

    SamplerBase::SampleLocType locs;
    cone_axis_locations(locs);
    count[1] = num_output_points();
    auto xyz = grp.var("points");
    xyz.put(locs[0].data(), start, count);

    if (m_output_cone_points) {
        std::vector<size_t> cone_start{nt, 0, 0};
        std::vector<size_t> cone_count{1, 0, AMREX_SPACEDIM};

        SamplerBase::SampleLocType conelocs;
        sampling_locations(conelocs);
        cone_count[1] = num_points();
        auto cone_xyz = grp.var("conepoints");
        cone_xyz.put(conelocs[0].data(), cone_start, cone_count);
    }
}

bool RadarSampler::output_netcdf_field(
    const std::vector<double>& /*output_buffer*/,
    ncutils::NCGroup& grp,
    const size_t nt)
{
    // Note: output_buffer is entire buffer...all samplers all
    // variables for this timestep

    m_los_velocity_interp.resize(m_los_velocity.size());

    for (int k = 0; k < m_ntotal; k++) {
        double mrat = static_cast<double>(k) / static_cast<double>(m_ntotal);
        for (int i = 0; i < m_npts; ++i) {
            int pi = i + k * m_npts;
            m_los_velocity_interp[pi] = mrat * m_los_velocity[pi] +
                                        (1.0 - mrat) * m_los_velocity_prior[pi];
        }
    }

    std::vector<size_t> start{nt, 0};
    std::vector<size_t> count{1, 0};
    start[1] = 0;
    count[1] = 0;

    auto var = grp.var("los_velocity");
    count[1] = num_output_points();
    var.put(m_los_velocity_interp.data(), start, count);

    return true;
}

#else

void RadarSampler::define_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
void RadarSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
void RadarSampler::output_netcdf_data(
    const ncutils::NCGroup& /*unused*/, const size_t /*unused*/) const
{}

bool RadarSampler::output_netcdf_field(
    const std::vector<double>& /*unused*/,
    ncutils::NCGroup& /*unused*/,
    const size_t /*unused*/)
{
    return true;
}

#endif

} // namespace amr_wind::sampling
