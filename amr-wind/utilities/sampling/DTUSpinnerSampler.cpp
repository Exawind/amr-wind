#include "amr-wind/utilities/sampling/DTUSpinnerSampler.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "AMReX_ParmParse.H"

#ifdef AMR_WIND_USE_OPENFAST
#include "amr-wind/wind_energy/actuator/Actuator.H"
#include "amr-wind/wind_energy/actuator/turbine/fast/TurbineFast.H"
#include "amr-wind/wind_energy/actuator/turbine/fast/turbine_fast_ops.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"
#endif

namespace amr_wind {
namespace sampling {

DTUSpinnerSampler::DTUSpinnerSampler(const CFDSim& sim) : LidarSampler(sim) {}

void DTUSpinnerSampler::initialize(const std::string& key)
{
    // Initialize the sampling time to be the same as simulation time
    m_time_sampling = m_sim.time().current_time();

    /*
     * Use this as a guide to implement the reading of the input
     * variables for the scanning pattern
     */

    // Read in new inputs specific to this class
    amrex::ParmParse pp(key);

    // Inner prism initial theta
    pp.get("inner_prism_theta0", m_InnerPrism.theta0);

    // Inner prism rotation rate
    pp.get("inner_prism_rotrate", m_InnerPrism.rot);

    // Inner prism azimuth angle
    pp.get("inner_prism_azimuth", m_InnerPrism.azimuth);

    // Outer prism initial theta
    pp.get("outer_prism_theta0", m_OuterPrism.theta0);

    // Outer prism rotation rate
    pp.get("outer_prism_rotrate", m_OuterPrism.rot);

    // Outer prism azimuth angle
    pp.get("outer_prism_azimuth", m_OuterPrism.azimuth);

    // This is the center of the lidar scan (x, y, z) [m]
    pp.getarr("lidar_center", m_lidar_center);
    AMREX_ALWAYS_ASSERT(
        static_cast<int>(m_lidar_center.size()) == AMREX_SPACEDIM);

    // Scan time
    pp.get("scan_time", m_scan_time);

    // Number of samples per scan
    pp.get("num_samples", m_num_samples);

    // Beam length
    pp.get("beam_length", m_beam_length);

    // Beam points
    pp.query("beam_points", m_beam_points);

    // Hub yaw logical flag
    pp.query("fixed_yaw", m_fixed_yaw);

    // Hub roll logical flag
    pp.query("fixed_roll", m_fixed_roll);

    // Hub tilt logical flag
    pp.query("fixed_tilt", m_fixed_tilt);

    // Spinner mode, hub or fixed
    pp.query("mode", m_spinner_mode);

    // For hub-mounted, this is turbine label
    pp.query("turbine", m_turbine_label);

    // Print hub fast turbine values to debug
    pp.query("hub_debug", m_hub_debug);

#ifdef AMR_WIND_USE_OPENFAST

    if (m_spinner_mode == "hub") {
        amrex::Print() << "Spinner Lidar will be attached to OpenFAST actuator"
                       << std::endl;
    }

#else

    AMREX_ALWAYS_ASSERT(m_spinner_mode == "fixed");

#endif // AMR_WIND_USE_OPENFAST

    update_sampling_locations();
    check_bounds();
}

vs::Vector DTUSpinnerSampler::adjust_lidar_pattern(
    vs::Vector beamPt, double yaw, double pitch, double roll)
{

    const vs::Vector angles(roll, pitch, yaw);
    auto beamPt_transform = sampling_utils::rotation(angles, beamPt);
    return beamPt_transform;
}

vs::Vector DTUSpinnerSampler::generate_lidar_pattern(
    PrismParameters InnerPrism, PrismParameters OuterPrism, double time)
{
    vs::Vector axis(1, 0, 0);
    vs::Vector ground(0, 0, 1);

    double innerTheta = InnerPrism.theta0 + InnerPrism.rot * time * 360;
    double outerTheta = OuterPrism.theta0 + OuterPrism.rot * time * 360;

    // NOLINTBEGIN(readability-suspicious-call-argument)
    auto reflection_1 = sampling_utils::rotate_euler_vec(
        axis, innerTheta,
        sampling_utils::rotate_euler_vec(
            ground, -(InnerPrism.azimuth / 2 + 90), axis));

    auto reflection_2 = sampling_utils::rotate_euler_vec(
        axis, outerTheta,
        sampling_utils::rotate_euler_vec(ground, OuterPrism.azimuth / 2, axis));
    // NOLINTEND(readability-suspicious-call-argument)

    return sampling_utils::reflect(
        reflection_2, sampling_utils::reflect(reflection_1, axis));
}

//

void DTUSpinnerSampler::sampling_locations(SampleLocType& locs) const
{

    // The total number of points at this time step
    long n_samples = m_beam_points * m_ntotal;

    // Resize to number of points in line times number of sampling times
    if (locs.size() < n_samples) {
        locs.resize(n_samples);
    }

    const amrex::Real ndiv = amrex::max(m_beam_points - 1, 1);
    amrex::Array<amrex::Real, AMREX_SPACEDIM> dx;

    // Loop per subsampling
    for (int k = 0; k < m_ntotal; ++k) {

        int offset = k * AMREX_SPACEDIM;

        // Loop per spatial dimension
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            dx[d] = (m_end[d + offset] - m_start[d + offset]) / ndiv;
        }

        for (int i = 0; i < m_beam_points; ++i) {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                locs[i + k * m_beam_points][d] =
                    m_start[d + offset] + i * dx[d];
            }
        }
    }
}

#ifdef AMR_WIND_USE_OPENFAST
void DTUSpinnerSampler::bcast_turbine(double* turbine_pack, int root_proc)
{
    BL_PROFILE("amr-wind::Sampling::DTUSpinnerSampler::bcast_turbine");

    amrex::ParallelDescriptor::Bcast(
        turbine_pack, 18, root_proc, amrex::ParallelDescriptor::Communicator());

    for (int i = 0; i < 9; i++) {
        current_hub_orient[i] = turbine_pack[i];
        if (i < 3) {
            current_hub_abs_pos[i] = turbine_pack[i + 9];
            current_hub_rot_vel[i] = turbine_pack[i + 12];
            turbine_base_pos[i] = turbine_pack[i + 15];
        }
    }
}

void DTUSpinnerSampler::get_turbine_data(std::string turbine_label)
{

    BL_PROFILE("amr-wind::Sampling::DTUSpinnerSampler::get_turbine_data");

    // Use Physics Manager to get actuators
    const auto& all_actuators =
        m_sim.physics_manager().get<actuator::Actuator>();
    actuator::ActuatorModel& lidar_act =
        all_actuators.get_act_bylabel(turbine_label);

    if (m_sim.time().current_time() == m_sim.time().start_time()) {
        std::string actlabel = lidar_act.label();
        amrex::Print() << "Spinner Lidar Attached to actuator: " << actlabel
                       << std::endl;
    }

    // Test cast to Line or Disk
    actuator::ActModel<actuator::TurbineFast, actuator::ActSrcLine>* actline =
        dynamic_cast<
            actuator::ActModel<actuator::TurbineFast, actuator::ActSrcLine>*>(
            &lidar_act);

    actuator::ActModel<actuator::TurbineFast, actuator::ActSrcDisk>* actdisk =
        dynamic_cast<
            actuator::ActModel<actuator::TurbineFast, actuator::ActSrcDisk>*>(
            &lidar_act);

    bool testline{(actline != nullptr)};
    bool testdisk{(actdisk != nullptr)};

    // Read and broadcast data
    if (testline) {
        auto& actdata = actline->meta().fast_data;
        const auto& info = actline->info();

        // Create buffer object
        double turbine_pack[18] = {};

        // Pack, broadcast, then unpack
        for (int i = 0; i < 9; i++) {
            turbine_pack[i] = actdata.hub_orient[i];
            if (i < 3) {
                turbine_pack[i + 9] = actdata.hub_abs_pos[i];
                turbine_pack[i + 12] = actdata.hub_rot_vel[i];
                turbine_pack[i + 15] = actdata.base_pos[i];
            }
        }

        bcast_turbine(turbine_pack, info.root_proc);
    } else if (testdisk) {
        auto& actdata = actdisk->meta().fast_data;
        const auto& info = actdisk->info();

        // Create buffer object
        double turbine_pack[18] = {};

        // Pack, broadcast, then unpack
        for (int i = 0; i < 9; i++) {
            turbine_pack[i] = actdata.hub_orient[i];
            if (i < 3) {
                turbine_pack[i + 9] = actdata.hub_abs_pos[i];
                turbine_pack[i + 12] = actdata.hub_rot_vel[i];
                turbine_pack[i + 15] = actdata.base_pos[i];
            }
        }

        bcast_turbine(turbine_pack, info.root_proc);
    } else {
        amrex::Abort("DTUSpinnerSampler: Problem finding actuator");
    }
}
#endif // AMR_WIND_USE_OPENFAST

void DTUSpinnerSampler::update_sampling_locations()
{
    BL_PROFILE(
        "amr-wind::Sampling::DTUSpinnerSampler::update_sampling_locations");

#ifdef AMR_WIND_USE_OPENFAST
    if (m_spinner_mode == "hub") {
        get_turbine_data(m_turbine_label);

        m_hub_location = vs::Vector(
            current_hub_abs_pos[0] + turbine_base_pos[0],
            current_hub_abs_pos[1] + turbine_base_pos[1],
            current_hub_abs_pos[2] + turbine_base_pos[2]);

        // TODO: Do we need an offset from the hub location
        // to lidar start along shaft axis? Same for static angle misalignment?

        m_lidar_center[0] = m_hub_location[0];
        m_lidar_center[1] = m_hub_location[1];
        m_lidar_center[2] = m_hub_location[2];

        m_hub_tilt = std::atan2(
            -current_hub_orient[6], std::sqrt(
                                        std::pow(current_hub_orient[7], 2.0) +
                                        std::pow(current_hub_orient[8], 2.0)));
        m_hub_roll = std::atan2(current_hub_orient[7], current_hub_orient[8]);
        m_hub_yaw = std::atan2(current_hub_orient[3], current_hub_orient[0]);
    }
#endif

    // Sampling called in post_advance so time should be new_time
    // Not current_time()
    amrex::Real time = m_sim.time().new_time();
    amrex::Real start_time = m_sim.time().start_time();
    amrex::Real dt_sim = m_sim.time().deltaT();
    const amrex::Real dt_s = m_scan_time / m_num_samples;
    amrex::Real start_diff = std::abs(time - start_time);

    // Initialize the sampling time to the first time in the simulation
    if (start_diff < 1e-10 && m_update_count == 0) {
        m_time_sampling = time;
        m_hub_location_init = m_hub_location;
    }

    amrex::Real ts_diff = time - m_time_sampling;

    // Correction for time mismatch
    int time_corr = (ts_diff > dt_s) ? int(ts_diff / dt_s) : 0;

    m_ns = int(dt_sim / dt_s) + time_corr;
    m_ntotal = int(std::ceil(dt_sim / dt_s));

    int n_size = AMREX_SPACEDIM * (m_ns);
    int n_totalsize = AMREX_SPACEDIM * (m_ntotal);

    if (m_hub_debug == true) {
        amrex::Print() << "ts_diff: " << ts_diff << "\t"
                       << "m_ns: " << m_ns << "\t"
                       << "Spin Time: " << m_time_sampling << "\t"
                       << "AMR Time: " << time << "\t"
                       << "n_size: " << n_size << std::endl;
    }

    // Resize these variables so they can store all the locations
    if (m_start.size() < n_totalsize) {
        m_start.resize(n_totalsize);
    }
    if (m_end.size() < n_totalsize) {
        m_end.resize(n_totalsize);
    }

#ifdef AMR_WIND_USE_OPENFAST
    if (m_hub_debug == true) {
        amrex::Print() << "Turbine Hub Pos: " << current_hub_abs_pos[0] << " "
                       << current_hub_abs_pos[1] << " "
                       << current_hub_abs_pos[2] << " " << std::endl;
        amrex::Print() << "Turbine Rot Vel: " << current_hub_rot_vel[0] << " "
                       << current_hub_rot_vel[1] << " "
                       << current_hub_rot_vel[2] << " " << std::endl;
        amrex::Print() << "Turbine Orient: " << current_hub_orient[0] << " "
                       << current_hub_orient[1] << " " << current_hub_orient[2]
                       << " " << std::endl;
        amrex::Print() << "Turbine Orient: " << current_hub_orient[3] << " "
                       << current_hub_orient[4] << " " << current_hub_orient[5]
                       << " " << std::endl;
        amrex::Print() << "Turbine Orient: " << current_hub_orient[6] << " "
                       << current_hub_orient[7] << " " << current_hub_orient[8]
                       << " " << std::endl;
        amrex::Print() << "Yaw:Tilt:Roll: " << m_hub_yaw << " " << m_hub_tilt
                       << " " << m_hub_roll << std::endl;
        amrex::Print() << "Last Yaw:Tilt:Roll: " << m_last_hub_yaw << " "
                       << m_last_hub_tilt << " " << m_last_hub_roll
                       << std::endl;
    }
#endif

    if (m_update_count > 0) {
        // Loop per subsampling
        for (int k = 0; k < m_ntotal; k++) {

            int offset = k * AMREX_SPACEDIM;

            if (k < m_ns) {
                amrex::Real step = (start_diff > 1e-10) ? 1.0 * k : 0.0;
                amrex::Real srat = (start_diff > 1e-10) ? step / m_ns : 0.0;

                // Unit vector in the direction of the beam
                auto beam_vector = generate_lidar_pattern(
                    m_InnerPrism, m_OuterPrism, m_time_sampling);

                // Use fancy trick to fix interp at 180 to -180 transition of
                // atan2
                amrex::Real step_hub_yaw =
                    m_last_hub_yaw +
                    (std::fmod(
                         std::fmod(m_hub_yaw - m_last_hub_yaw, twopi) + threepi,
                         twopi) -
                     pi) *
                        srat;
                amrex::Real step_hub_tilt =
                    m_last_hub_tilt +
                    (std::fmod(
                         std::fmod(m_hub_tilt - m_last_hub_tilt, twopi) +
                             threepi,
                         twopi) -
                     pi) *
                        srat;
                amrex::Real step_hub_roll =
                    m_last_hub_roll +
                    (std::fmod(
                         std::fmod(m_hub_roll - m_last_hub_roll, twopi) +
                             threepi,
                         twopi) -
                     pi) *
                        srat;

                // Rotate beam unit vector
                beam_vector = adjust_lidar_pattern(
                    beam_vector, m_fixed_yaw + step_hub_yaw * radtodeg,
                    m_fixed_tilt + step_hub_tilt * radtodeg,
                    m_fixed_roll + step_hub_roll * radtodeg);

                // Interpolate lidar center
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    step_lidar_center[d] =
                        (step / m_ns) *
                            (m_lidar_center[d] - m_last_lidar_center[d]) +
                        m_last_lidar_center[d];
                }

                // Beam start and end points
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    m_start[d + offset] = step_lidar_center[d];
                    m_end[d + offset] =
                        step_lidar_center[d] + beam_vector[d] * m_beam_length;
                }

                m_time_sampling += dt_s;
            } else {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    m_start[d + offset] = -99999.99;
                    m_end[d + offset] = -99999.99;
                }
            }
        }
    }

    m_last_hub_yaw = m_hub_yaw;
    m_last_hub_tilt = m_hub_tilt;
    m_last_hub_roll = m_hub_roll;
    m_last_lidar_center = m_lidar_center;

    m_update_count = m_update_count + 1;
}

#ifdef AMR_WIND_USE_NETCDF

void DTUSpinnerSampler::define_netcdf_metadata(
    const ncutils::NCGroup& grp) const
{
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("start", m_start);
    grp.put_attr("end", m_end);
    grp.put_attr("time_table", m_time_table);
    grp.put_attr("azimuth_table", m_azimuth_table);
    grp.put_attr("elevation_table", m_elevation_table);

    std::vector<double> scan_time{m_scan_time};
    std::vector<double> num_samples{m_num_samples};
    std::vector<double> beam_length{m_beam_length};
    std::vector<int> beam_points{m_beam_points};
    grp.put_attr("scan_time", scan_time);
    grp.put_attr("num_samples", num_samples);
    grp.put_attr("beam_length", beam_length);
    grp.put_attr("beam_points", beam_points);

    std::vector<double> fixed_tilt{m_fixed_tilt};
    std::vector<double> fixed_yaw{m_fixed_yaw};
    std::vector<double> fixed_roll{m_fixed_roll};
    grp.put_attr("fixed_tilt", fixed_tilt);
    grp.put_attr("fixed_yaw", fixed_yaw);
    grp.put_attr("fixed_roll", fixed_roll);

    std::vector<double> innerprism_theta0{m_InnerPrism.theta0};
    std::vector<double> innerprism_rot{m_InnerPrism.rot};
    std::vector<double> innerprism_azimuth{m_InnerPrism.azimuth};
    std::vector<double> outerprism_theta0{m_OuterPrism.theta0};
    std::vector<double> outerprism_rot{m_OuterPrism.rot};
    std::vector<double> outerprism_azimuth{m_OuterPrism.azimuth};
    grp.put_attr("innerprism_theta0", innerprism_theta0);
    grp.put_attr("innerprism_rot", innerprism_rot);
    grp.put_attr("innerprism_azimuth", innerprism_azimuth);
    grp.put_attr("outerprism_theta0", outerprism_theta0);
    grp.put_attr("outerprism_rot", outerprism_rot);
    grp.put_attr("outerprism_azimuth", outerprism_azimuth);

    grp.put_attr("spinner_mode", m_spinner_mode);
    grp.put_attr("turbine", m_turbine_label);

    grp.def_dim("ndim", AMREX_SPACEDIM);
    grp.def_var("points", NC_DOUBLE, {"num_time_steps", "num_points", "ndim"});

    grp.def_var("points_x", NC_DOUBLE, {"num_time_steps", "num_points"});
    grp.def_var("points_y", NC_DOUBLE, {"num_time_steps", "num_points"});
    grp.def_var("points_z", NC_DOUBLE, {"num_time_steps", "num_points"});

    grp.def_dim("nang", 3);
    grp.def_var("rotor_angles_rad", NC_DOUBLE, {"num_time_steps", "nang"});

    grp.def_var("rotor_hub_pos", NC_DOUBLE, {"num_time_steps", "ndim"});

    std::vector<double> fillnan{std::nan("1")};
    grp.var("rotor_hub_pos").put_attr("_FillValue", fillnan);
    grp.var("rotor_angles_rad").put_attr("_FillValue", fillnan);
    grp.var("points").put_attr("_FillValue", fillnan);
}

void DTUSpinnerSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}

void DTUSpinnerSampler::output_netcdf_data(
    const ncutils::NCGroup& grp, const size_t nt) const
{
    std::vector<size_t> start{nt, 0, 0};
    std::vector<size_t> count{1, 0, AMREX_SPACEDIM};
    std::vector<size_t> starti{nt, 0};
    std::vector<size_t> counti{1, 0};
    SamplerBase::SampleLocType locs;
    sampling_locations(locs);

    auto xyz = grp.var("points");
    auto xp = grp.var("points_x");
    auto yp = grp.var("points_y");
    auto zp = grp.var("points_z");
    count[1] = num_points();
    counti[1] = num_points();

    xyz.put(&locs[0][0], start, count);

    int n_samples = m_beam_points * m_ntotal;

    double xlocs[n_samples];
    double ylocs[n_samples];
    double zlocs[n_samples];

    // Loop per subsampling
    for (int k = 0; k < m_ntotal; ++k) {

        int offset = k * AMREX_SPACEDIM;

        for (int i = 0; i < m_beam_points; ++i) {
            xlocs[i + k * m_beam_points] = locs[i + k * m_beam_points][0];
            ylocs[i + k * m_beam_points] = locs[i + k * m_beam_points][1];
            zlocs[i + k * m_beam_points] = locs[i + k * m_beam_points][2];
        }
    }

    xp.put(&xlocs[0], starti, counti);
    yp.put(&ylocs[0], starti, counti);
    zp.put(&zlocs[0], starti, counti);

    auto angs = grp.var("rotor_angles_rad");
    std::vector<size_t> acount{1, 3};
    double rangs[3] = {m_hub_tilt, m_hub_roll, m_hub_yaw};
    angs.put(rangs, start, acount);

    auto hpos = grp.var("rotor_hub_pos");
    std::vector<size_t> pcount{1, AMREX_SPACEDIM};
    double rhpos[3] = {m_lidar_center[0], m_lidar_center[1], m_lidar_center[2]};
    hpos.put(rhpos, start, pcount);
}

#else

void DTUSpinnerSampler::define_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}

void DTUSpinnerSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}

void DTUSpinnerSampler::output_netcdf_data(
    const ncutils::NCGroup& /*unused*/, const size_t /*unused*/) const
{}

#endif

} // namespace sampling

template struct ::amr_wind::sampling::SamplerBase::Register<
    ::amr_wind::sampling::DTUSpinnerSampler>;

} // namespace amr_wind
