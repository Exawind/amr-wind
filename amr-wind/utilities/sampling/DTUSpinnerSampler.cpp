#include "amr-wind/utilities/sampling/DTUSpinnerSampler.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "amr-wind/wind_energy/actuator/Actuator.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"
#include "amr-wind/wind_energy/actuator/turbine/fast/TurbineFast.cpp"

#include "AMReX_ParmParse.H"

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

    // Turbine yaw angle
    pp.get("turbine_yaw_angle", m_turbine_yaw_angle);

    // Hub yaw logical flag
    pp.query("hub_yaw", m_hub_yaw);

    // Hub roll logical flag
    pp.query("hub_roll", m_hub_roll);

    // Hub tilt logical flag
    pp.query("hub_tilt", m_hub_tilt);

    // Spinner mode, hub or fixed
    pp.query("mode", m_spinner_mode);

    // For hub-mounted, this is turbine label
    pp.query("turbine", m_turbine_label);

    // Print hub fast turbine values to debug
    pp.query("hub_debug", m_hub_debug);

    update_sampling_locations();

    check_bounds();
}

// functions
auto reflect(vs::Vector line, vs::Vector vec)
{

    vs::Tensor ref(
        1 - 2 * line.x() * line.x(), -2 * line.x() * line.y(),
        -2 * line.x() * line.z(), -2 * line.y() * line.x(),
        1 - 2 * line.y() * line.y(), -2 * line.y() * line.z(),
        -2 * line.z() * line.x(), -2 * line.z() * line.y(),
        1 - 2 * line.z() * line.z());

    return vec & ref;
}

auto rotate_euler_vec(vs::Vector axis, double angle, vs::Vector vec)
{

    axis.normalize();

    const auto RotMat = vs::quaternion(axis, angle);

    return vec & RotMat;
}

vs::Vector rotation(const vs::Vector& angles, const vs::Vector& data)
{
    const vs::Tensor rotMatrix =
        vs::xrot(angles.x()) & vs::yrot(angles.y()) & vs::zrot(angles.z());
    return data & rotMatrix;
}

vs::Vector
adjust_lidar_pattern(vs::Vector beamPt, double yaw, double pitch, double roll)
{

    const vs::Vector angles(roll, pitch, yaw);

    auto beamPt_transform = rotation(angles, beamPt);

    return beamPt_transform;
}

vs::Vector generate_lidar_pattern(
    PrismParameters InnerPrism, PrismParameters OuterPrism, double time)
{
    vs::Vector axis(1, 0, 0);

    const vs::Vector ground(0, 0, 1);

    const double innerTheta = InnerPrism.theta0 + InnerPrism.rot * time * 360;
    const double outerTheta = OuterPrism.theta0 + OuterPrism.rot * time * 360;

    const auto reflection_1 = rotate_euler_vec(
        axis, innerTheta,
        rotate_euler_vec(ground, -(InnerPrism.azimuth / 2 + 90), axis));

    const auto reflection_2 = rotate_euler_vec(
        axis, outerTheta,
        rotate_euler_vec(ground, OuterPrism.azimuth / 2, axis));

    return reflect(reflection_2, reflect(reflection_1, axis));
}

//

void DTUSpinnerSampler::sampling_locations(SampleLocType& locs) const
{

    // The total number of points at this time step
    int n_samples = m_beam_points * m_ns;

    // Resize to number of points in line times number of sampling times
    if (locs.size() < n_samples) {
        locs.resize(n_samples);
    }

    const amrex::Real ndiv = amrex::max(m_beam_points - 1, 1);
    amrex::Array<amrex::Real, AMREX_SPACEDIM> dx;

    // Loop per subsampling
    for (int k = 0; k < m_ns; ++k) {

        int offset = k * AMREX_SPACEDIM;

        // Loop per spacial dimension
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

void DTUSpinnerSampler::bcast_turbine(
    const exw_fast::FastTurbine& actdata, int root_proc)
{
    BL_PROFILE("amr-wind::Sampling::DTUSpinnerSampler::bcast_turbine");

    // Create buffer object
    amrex::Real turbine_pack[18] = {};

    // Pack, broadcast, then unpack
    for (int i = 0; i < 9; i++) {
        turbine_pack[i] = actdata.hub_orient[i];
        if (i < 3) {
            turbine_pack[i + 9] = actdata.hub_abs_pos[i];
            turbine_pack[i + 12] = actdata.hub_rot_vel[i];
            turbine_pack[i + 15] = actdata.base_pos[i];
        }
    }

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
        nullptr;
    actline = dynamic_cast<
        actuator::ActModel<actuator::TurbineFast, actuator::ActSrcLine>*>(
        &lidar_act);

    actuator::ActModel<actuator::TurbineFast, actuator::ActSrcDisk>* actdisk =
        nullptr;
    actdisk = dynamic_cast<
        actuator::ActModel<actuator::TurbineFast, actuator::ActSrcDisk>*>(
        &lidar_act);

    bool testline{(actline != nullptr)};
    bool testdisk{(actdisk != nullptr)};

    // Read and broadcast data
    if (testline) {
        auto& actdata = actline->meta().fast_data;
        const auto& info = actline->info();
        bcast_turbine(actdata, info.root_proc);
    } else if (testdisk) {
        auto& actdata = actdisk->meta().fast_data;
        const auto& info = actdisk->info();
        bcast_turbine(actdata, info.root_proc);
    } else {
        amrex::Abort("DTUSpinnerSampler: Problem finding actuator");
    }
}

void DTUSpinnerSampler::update_sampling_locations()
{
    BL_PROFILE(
        "amr-wind::Sampling::DTUSpinnerSampler::update_sampling_locations");

    get_turbine_data(m_turbine_label);

    if (m_spinner_mode == "hub") {
        m_hub_location = vs::Vector(
            current_hub_abs_pos[0] + turbine_base_pos[0],
            current_hub_abs_pos[1] + turbine_base_pos[1],
            current_hub_abs_pos[2] + turbine_base_pos[2]);

        // TODO: Do we need an offset from the hub location
        // to lidar start along shaft axis? Same for static angle misalignment?

        m_lidar_center[0] = m_hub_location[0];
        m_lidar_center[1] = m_hub_location[1];
        m_lidar_center[2] = m_hub_location[2];

        m_hub_tilt = -std::atan2(
                         -current_hub_orient[6],
                         std::sqrt(
                             std::pow(current_hub_orient[7], 2.0) +
                             std::pow(current_hub_orient[8], 2.0))) *
                     180.0 / M_PI;
        m_hub_roll = std::atan2(current_hub_orient[7], current_hub_orient[8]) *
                     180.0 / M_PI;
        m_hub_yaw = std::atan2(current_hub_orient[3], current_hub_orient[0]) *
                    180.0 / M_PI;
    }

    amrex::Real time = m_sim.time().current_time();
    amrex::Real start_time = m_sim.time().start_time();
    amrex::Real dt_sim = m_sim.time().deltaT();

    // Initialize the sampling time to the first time in the simulation
    if (time == start_time) {
        m_time_sampling = time;
        m_hub_location_init = m_hub_location;
    }

    // Timestep of spinner lidar
    const amrex::Real dt_s = m_scan_time / m_num_samples;

    // Integer number of spinner points per timestep
    m_ns = int(dt_sim / dt_s);
    int n_size = AMREX_SPACEDIM * m_ns;

    // Resize these variables so they can store all the locations
    if (m_start.size() < n_size) {
        m_start.resize(n_size);
    }
    if (m_end.size() < n_size) {
        m_end.resize(n_size);
    }

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

    // Loop per subsampling
    for (int k = 0; k < m_ns; ++k) {

        int offset = k * AMREX_SPACEDIM;

        m_time_sampling += dt_s;

        float step = static_cast<float>(k) + 1.0;

        // Unit vector in the direction of the beam
        auto beam_vector =
            generate_lidar_pattern(m_InnerPrism, m_OuterPrism, m_time_sampling);

        // Interpolate angles over sim time
        amrex::Real step_hub_yaw =
            (step / m_ns) * (m_hub_yaw - m_last_hub_yaw) + m_last_hub_yaw;
        amrex::Real step_hub_tilt =
            (step / m_ns) * (m_hub_tilt - m_last_hub_tilt) + m_last_hub_tilt;
        amrex::Real step_hub_roll =
            (step / m_ns) * (m_hub_roll - m_last_hub_roll) + m_last_hub_roll;

        // Rotate beam unit vector
        beam_vector = adjust_lidar_pattern(
            beam_vector, step_hub_yaw, step_hub_tilt, step_hub_roll);

        // Interpolate lidar center
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            step_lidar_center[d] =
                (step / m_ns) * (m_lidar_center[d] - m_last_lidar_center[d]) +
                m_last_lidar_center[d];
        }

        // Beam start and end points
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            m_start[d + offset] = step_lidar_center[d];
            m_end[d + offset] =
                step_lidar_center[d] + beam_vector[d] * m_beam_length;
        }
    }

    m_last_hub_yaw = m_hub_yaw;
    m_last_hub_tilt = m_hub_tilt;
    m_last_hub_roll = m_hub_roll;
    m_last_lidar_center = m_lidar_center;
}

#ifdef AMR_WIND_USE_NETCDF

bool DTUSpinnerSampler::output_netcdf_field(
    double* /*unused*/, ncutils::NCVar& /*unused*/)
{
    return true;
}

void DTUSpinnerSampler::define_netcdf_metadata(
    const ncutils::NCGroup& grp) const
{
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("start", m_start);
    grp.put_attr("end", m_end);
    grp.put_attr("time_table", m_time_table);
    grp.put_attr("azimuth_table", m_azimuth_table);
    grp.put_attr("elevation_table", m_elevation_table);
    grp.def_var("points", NC_DOUBLE, {"num_time_steps", "num_points", "ndim"});
}

void DTUSpinnerSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}

void DTUSpinnerSampler::output_netcdf_data(
    const ncutils::NCGroup& grp, const size_t nt) const
{
    // Write the coordinates every time
    std::vector<size_t> start{nt, 0, 0};
    std::vector<size_t> count{1, 0, AMREX_SPACEDIM};
    SamplerBase::SampleLocType locs;
    sampling_locations(locs);
    auto xyz = grp.var("points");
    count[1] = num_points();
    xyz.put(&locs[0][0], start, count);
}

#else

bool DTUSpinnerSampler::output_netcdf_field(
    double* /*unused*/, ncutils::NCVar& /*unused*/)
{
    return true;
}

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
