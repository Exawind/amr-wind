#include "amr-wind/utilities/sampling/RadarSampler.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tensor_ops.H"

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
    pp.getarr("endpoint", m_end);
    //AMREX_ALWAYS_ASSERT(static_cast<int>(m_origin.size()) == AMREX_SPACEDIM);

    // Sampling Frequency 
    pp.get("sampling_frequency", m_sample_freq);

    // Sampling Frequency 
    pp.get("radar_sampling_frequency", m_radar_sample_freq);

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

    // Unit vector axis of the beam [-]
    pp.getarr("axis", m_axis);

    // Unit vertical axis [-]
    pp.getarr("vertical_unit_dir", m_vertical);

    // Sweep angle of the beam [deg]
    pp.get("sweep_angle", m_sweep_angle);

    // Reset time of the beam at elevation change [s]
    pp.get("reset_time", m_reset_time);

    // Elevation angles [deg]
    pp.getarr("elevation_angles", m_elevation_angles);

    // Set initial periodic time to half-way point of forward sweep
    m_periodic_time = 0.5*m_sweep_angle / m_angular_speed;
    m_current_phase = determine_operation_phase();

    m_cone_size = num_points_cone();

    RadarSampler::new_cone();

    amrex::Print() << "Passed originate cone" << std::endl;

    RadarSampler::update_sampling_locations();

    amrex::Print() << "Passed update sampling" << std::endl;

    check_bounds();

    amrex::Print() << "Passed check bounds" << std::endl;
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
  m_periodic_time = m_radar_time - std::floor(m_radar_time / total_sweep_time()) * total_sweep_time();
  return m_periodic_time;
}

// Single direction sweep count
int RadarSampler::sweep_count() const
{
  const double sweep_time = m_sweep_angle / m_angular_speed + m_reset_time;
  return std::floor(m_radar_time / sweep_time);
}

RadarSampler::phase 
RadarSampler::determine_operation_phase() const
{
  const double phase_time = m_sweep_angle / m_angular_speed;

  if (m_periodic_time < phase_time) {
    return phase::FORWARD;
  } else if (m_periodic_time < phase_time + m_reset_time) {
    return phase::FORWARD_PAUSE;
  } else if (m_periodic_time < 2 * phase_time + m_reset_time) {
    return phase::REVERSE;
  } else {
    return phase::REVERSE_PAUSE;
  }
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

void RadarSampler::update_sampling_locations()
{

    //amrex::Print() << "RADAR UPDATE SAMPLING LOCATIONS" << std::endl;
    amrex::Real time = m_sim.time().current_time();
    amrex::Real start_time = m_sim.time().start_time();
    amrex::Real dt_sim = m_sim.time().deltaT();
    amrex::Real dt_sample = 1.0/m_sample_freq;
    amrex::Real ts_diff = time - m_radar_time;

    if (time == start_time) {
        m_radar_time = time;
    }
    
    // Correction for time mismatch
    int time_corr = (ts_diff > dt_sample) ? int(ts_diff / dt_sample) : 0;

    m_ns = int(dt_sim / dt_sample) + time_corr;
    m_ntotal = int(std::ceil(dt_sim / dt_sample));

    int n_allcones = m_ntotal*num_points();
    current_cones.resize(n_allcones);

    amrex::Print() << "ts_diff: " << ts_diff << "\t" 
                   << "m_ns: " << m_ns << "\t"
                   << "m_ntotal: " << m_ntotal << "\t"
                   << "Radar Time: " << m_radar_time << "\t"
                   << "AMR Time: " << time << std::endl;

    // Loop for oversampling
    for (int k = 0; k < m_ntotal; k++) {
        if (k < m_ns) {
            double per_time = periodic_time();
            double sweep_angle = determine_current_sweep_angle();
            double elevation_angle = m_elevation_angles.at(sweep_count() % m_elevation_angles.size());

            amrex::Print() << "-------------------------" << std::endl
                        << "Total Sweep: " << total_sweep_time() << "\t"
                        << "Periodic Time: " << per_time << "\t" 
                        << "Radar Time: " << m_radar_time << "\t"
                        << "S Angle: " << sweep_angle << "\t"
                        << "E Angle: " << elevation_angle << "\t"
                        << "Sim Time: " << time << std::endl
                        << "-------------------------"  << std::endl;

            vs::Vector vertical_ref(m_vertical[0],m_vertical[1],m_vertical[2]);
            vs::Vector radar_ref(m_axis[0],m_axis[1],m_axis[2]);
            vs::Vector swept_axis(sampling_utils::rotate_euler_vector(vertical_ref,sweep_angle,radar_ref));
            vs::Vector normal_ref(-vertical_ref ^ swept_axis);
            //vs::Vector rotation_axis(sampling_utils::rotate_euler_vector(normal_ref,elevation_angle,vertical_ref));

            // Add rotated cone to current cones
            for (int i = 0; i < m_cone_size; ++i) {
                vs::Vector temp_point( initial_cone[i][0] - m_start[0], 
                                       initial_cone[i][1] - m_start[1], 
                                       initial_cone[i][2] - m_start[2]);
                vs::Vector swept_point(sampling_utils::rotate_euler_vector(vertical_ref,sweep_angle,temp_point));
                vs::Vector rotated_point(sampling_utils::rotate_euler_vector(normal_ref,elevation_angle,swept_point));
                int point_index = i + k*m_cone_size;
                current_cones[point_index][0] = rotated_point[0] + m_start[0];
                current_cones[point_index][1] = rotated_point[1] + m_start[1];
                current_cones[point_index][2] = rotated_point[2] + m_start[2];
            }

            m_radar_time += dt_sample;
        }else{
            // Fill with -99999.99
            for (int i = 0; i < m_cone_size; ++i) {
                int point_index = i + k*m_cone_size;
                current_cones[point_index][0] = -99999.99;
                current_cones[point_index][1] = -99999.99;
                current_cones[point_index][2] = -99999.99;
            }
        }
    }   
}

void RadarSampler::new_cone()
{
    vs::Vector canon_vector(0.0, 0.0, 1.0);
    vs::Vector origin_vector(0.0, 0.0, 0.0);
    vs::Vector init_axis(m_axis[0],m_axis[1],m_axis[2]);

    initial_cone.resize(num_points_cone());
    m_rays.resize(num_points_quad());
    m_weights.resize(num_points_quad());

    sampling_utils::spherical_cap_truncated_normal(m_cone_angle, ntheta, sampling_utils::NormalRule::HALFPOWER, m_rays, m_weights);

    int nquad = num_points_quad();

    for (int j = 0; j < nquad; ++j) {
      amrex::Print() << "Rays x: " << m_rays[j][0] << "\t" << "Rays y: " << m_rays[j][1] << "\t" << "Rays z: " << m_rays[j][2] << "\t" << std::endl;
    }

    amrex::Array<amrex::Real, AMREX_SPACEDIM> dx;

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        dx[d] = (canon_vector[d] - origin_vector[d]) / m_npts;
    }

    amrex::Print() << "-------------------------" << std::endl
                << "Beam Length: " << m_beam_length << "\t"
                << "nquad: " << nquad << "\t" 
                << "dx[2]: " << dx[2] << "\t"
                << "m_npts: " << m_npts << "\t"
                << "ntheta: " << ntheta << "\t"
                << "nphi: " << nphi << std::endl
                << "-------------------------"  << std::endl;

    for (int n = 0; n < m_npts; ++n) {
      for (int j = 0; j < nquad; ++j) {
        int pt_idx = j + n*nquad;
        const auto radius = dx[2]*n;
        initial_cone[pt_idx][0] = radius*m_rays[j][0];
        initial_cone[pt_idx][1] = radius*m_rays[j][1];
        initial_cone[pt_idx][2] = radius*m_rays[j][2];
      }
    }

    vs::Vector tc_axis = canon_vector ^ init_axis;
    amrex::Real tc_angle = std::acos(init_axis & canon_vector)*180.0/PI;
    
    for (int i = 0; i < num_points_cone(); ++i) {
        vs::Vector temp_loc( initial_cone[i][0], initial_cone[i][1], initial_cone[i][2]);
        vs::Vector new_rot(sampling_utils::rotate_euler_vector(tc_axis,tc_angle,temp_loc));
        initial_cone[i][0] = new_rot[0]*m_beam_length + m_start[0];
        initial_cone[i][1] = new_rot[1]*m_beam_length + m_start[1];
        initial_cone[i][2] = new_rot[2]*m_beam_length + m_start[2];
    }

}

void RadarSampler::origin_cone()
{

    vs::Vector canon_vector(0.0, 0.0, 1.0);
    vs::Vector origin_vector(0.0, 0.0, 0.0); 

    
    int full_size = m_npts + ntheta*(m_npts-1)*(nphi-1);
    initial_cone.resize(full_size);

    const amrex::Real nptsdiv = amrex::max(m_npts - 1, 1);
    const amrex::Real nthediv = amrex::max(ntheta, 1);
    const amrex::Real nphidiv = amrex::max(nphi - 1, 1);
    amrex::Array<amrex::Real, AMREX_SPACEDIM> dx;

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        dx[d] = (canon_vector[d] - origin_vector[d]) / nptsdiv;
    }

    // Assemble Radar axis unit beam
    for (int i = 0; i < m_npts; ++i) {
        //amrex::Print() << "i: " << i << std::endl;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            //locs[i][d] = origin_vector[d] + i * dx[d];
            initial_cone[i][d] = origin_vector[d] + i * dx[d];
        }
    }

    const amrex::Real dphi = m_cone_angle / nphidiv;
    const amrex::Real dtheta = 360.0 / nthediv;
    amrex::Real c_theta = 0.0;
    amrex::Real c_phi = 0.0;

    // Assemble all other unit beams
    for (int t = 0; t < ntheta; ++t) {
        c_theta = t * dtheta;
        for (int p = 0; p < nphi-1; ++p) {
            c_phi = static_cast<double>((p + 1) * dphi);
            for (int i = 0; i < m_npts-1; ++i) {
                int ci = m_npts + i + (m_npts-1)*(p + (nphi-1)*t);

                const vs::Vector temp_loc(origin_vector[0]+dx[0] + i * dx[0],
                                          origin_vector[1]+dx[1] + i * dx[1],
                                          origin_vector[2]+dx[2] + i * dx[2]);
                const vs::Vector temp_ang(0.0,c_phi,c_theta); 
                vs::Vector temp_output = sampling_utils::canon_rotator(temp_ang,temp_loc);

                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    initial_cone[ci][d] = temp_output[d];
                }
            }
        }
    }

    // Transform to scale and align Radar axis with start and end points
    //vs::Vector target_beam(m_end[0]-m_start[0],m_end[1]-m_start[1],m_end[2]-m_start[2]); 
    //amrex::Real target_beam_scale = std::sqrt(std::pow(target_beam[0],2.0)+
    //                    std::pow(target_beam[1],2.0)+
    //                    std::pow(target_beam[2],2.0));
    vs::Vector init_axis(m_axis[0],m_axis[1],m_axis[2]);
    vs::Vector tc_axis = canon_vector ^ init_axis;
    amrex::Real tc_angle = std::acos(init_axis & canon_vector)*180.0/PI;
    
    for (int i = 0; i < full_size; ++i) {
        vs::Vector temp_loc( initial_cone[i][0], initial_cone[i][1], initial_cone[i][2]);
        vs::Vector new_rot(sampling_utils::rotate_euler_vector(tc_axis,tc_angle,temp_loc));
        initial_cone[i][0] = new_rot[0]*m_beam_length + m_start[0];
        initial_cone[i][1] = new_rot[1]*m_beam_length + m_start[1];
        initial_cone[i][2] = new_rot[2]*m_beam_length + m_start[2];
    }

}

std::vector<double> RadarSampler::modify_sample_data(const std::vector<double>& sample_data)
{
  // sample_data enters this method for each sampled variable
  // there are m_ntotal steps (cones) based on device sampling rate
  AMREX_ALWAYS_ASSERT(static_cast<int>(sample_data.size()) == num_points());
  
  const int n_buf = sample_data.size();
  const int n_cones = m_ntotal;
  const int n_vars = n_buf/num_points();

  std::vector<double> mod_data(num_output_points());

  for (int iv = 0; iv < n_vars; ++iv) {
    for (int ic = 0; ic < n_cones; ++ic) {
      int c_idx = iv*n_cones + ic;
      int c_start = c_idx*num_points_cone();
      int c_end = (c_idx + 1)*num_points_cone();
      int a_start = c_idx*num_points_axis();

      // Send a single cone to be line averaged 
      std::vector<double> temp_vals(sample_data.begin() + c_start, sample_data.begin() + c_end);
      line_average(m_weights,temp_vals,mod_data,a_start);
    }
  }
  return mod_data;
}


void RadarSampler::line_average(
  const std::vector<double>& weights,
  const std::vector<double>& values,
  std::vector<double>& reduced,
  int offset)
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
      reduced[n+offset] = average / weight_sum;
    } else {
      reduced[n+offset] = 0;
    }
  }
}


void RadarSampler::sampling_locations(SampleLocType& locs) const
{
    locs.resize(num_points());

    for (int i = 0; i < num_points(); ++i) {
        locs[i][0] = current_cones[i][0];
        locs[i][1] = current_cones[i][1];
        locs[i][2] = current_cones[i][2];
    }
}

void RadarSampler::cone_axis_locations(SampleLocType& axis_locs) const
{
    axis_locs.resize(num_output_points());

    //locs = create_cone();
    for (int k = 0; k < m_ntotal; k++) {
      for (int i = 0; i < m_npts; ++i) {
        int pi = i + k*m_npts;
        int ci = i*num_points_quad() + k*num_points_cone();
        axis_locs[pi][0] = current_cones[ci][0];
        axis_locs[pi][1] = current_cones[ci][1]; 
        axis_locs[pi][2] = current_cones[ci][2];
      }
    }
}

#ifdef AMR_WIND_USE_NETCDF
void RadarSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
{
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("start", m_start);
    grp.put_attr("end", m_end);

    grp.def_dim("ndim", AMREX_SPACEDIM);
    grp.def_dim("num_points_cone", num_points());
    grp.def_var("points", NC_DOUBLE, {"num_time_steps", "num_points", "ndim"});
    grp.def_var("conepoints", NC_DOUBLE, {"num_time_steps", "num_points_cone", "ndim"});
}

void RadarSampler::populate_netcdf_metadata(const ncutils::NCGroup&) const {}

void RadarSampler::output_netcdf_data(
    const ncutils::NCGroup& grp, const size_t nt) const
{
    amrex::Print() << "Output Radar Sampler Netcdf Data" << std::endl;
    std::vector<size_t> start{nt, 0, 0};
    std::vector<size_t> count{1, 0, AMREX_SPACEDIM};
    
    std::vector<size_t> cone_start{nt, 0, 0};
    std::vector<size_t> cone_count{1, 0, AMREX_SPACEDIM};
    
    SamplerBase::SampleLocType locs;
    cone_axis_locations(locs);
    count[1] = num_output_points();
    auto xyz = grp.var("points");
    xyz.put(&locs[0][0], start, count);

    SamplerBase::SampleLocType conelocs;
    sampling_locations(conelocs);
    cone_count[1] = num_points();
    auto cone_xyz = grp.var("conepoints");
    cone_xyz.put(&conelocs[0][0], cone_start, cone_count);
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

#endif

} // namespace amr_wind::sampling
