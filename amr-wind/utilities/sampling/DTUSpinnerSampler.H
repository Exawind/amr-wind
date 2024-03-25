#ifndef DTUSPINNERSAMPLER_H
#define DTUSPINNERSAMPLER_H

#include "amr-wind/utilities/sampling/SamplingUtils.H"
#include "amr-wind/utilities/sampling/LidarSampler.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind::sampling {

struct PrismParameters
{
    amrex::Real theta0{90};
    amrex::Real rot{3.5};
    amrex::Real azimuth{15.2};
};

/** Sample data along a line that moves in time
 *  \ingroup sampling
 *
 *  Defines probes along a line defined by a time history from `m_azimuth_table`
 *  and `m_elevation_table` with a length of `m_length`. This determines the
 * `start` and `end` points at every time-step and divides the line into equal
 * segments as determined by `num_points`. The length of each segment is given
 * by \f[ m_length / (num_points - 1) \f]
 */
class DTUSpinnerSampler : public LidarSampler
{
public:
    // cppcheck-suppress duplInheritedMember
    static std::string identifier() { return "DTUSpinnerSampler"; }

    explicit DTUSpinnerSampler(const CFDSim& /*sim*/);

    /** Read user inputs and initialize the sampling object
     *
     *  \param key Prefix used to parse inputs from file
     */
    void initialize(const std::string& key) override;

    /** Update the position of the sampling locations in time
     *
     *
     */
    void sampling_locations(SampleLocType& /*locs*/) const override;

    static vs::Vector generate_lidar_pattern(
        PrismParameters InnerPrism, PrismParameters OuterPrism, double time);

    static vs::Vector adjust_lidar_pattern(
        vs::Vector beamPt, double yaw, double pitch, double roll);

#ifdef AMR_WIND_USE_OPENFAST

    void bcast_turbine(double* /*turbine_pack*/, int /*root_proc*/);
    void get_turbine_data(std::string turbine_label);

#endif

    /** Determine how the subsampling will be done
     *
     */
    void update_sampling_locations() override;

    void post_sample_actions() override{};

    long num_points() const override { return (m_beam_points * m_ntotal); }

    long num_output_points() const override
    {
        return (m_beam_points * m_ntotal);
    }

    void
    define_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;
    void
    populate_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;
    void output_netcdf_data(
        const ncutils::NCGroup& /*unused*/,
        const size_t /*unused*/) const override;

private:
    // Number of subsamples in this timestep
    int m_ns{1};

    // Max number of subsamples in all timesteps
    long m_ntotal{1};

    int m_update_count{0};
    const amrex::Real radtodeg{180.0 / M_PI};
    const amrex::Real pi{M_PI};
    const amrex::Real twopi{2.0 * M_PI};
    const amrex::Real threepi{3.0 * M_PI};
    amrex::Real m_time_sampling{0};
    PrismParameters m_InnerPrism;
    PrismParameters m_OuterPrism;
    amrex::Vector<amrex::Real> m_lidar_center{0.0, 0.0, 0.0};
    amrex::Real m_scan_time{2};
    amrex::Real m_num_samples{984};
    amrex::Real m_beam_length{270};
    int m_beam_points{432};
    amrex::Real m_fixed_yaw{0};
    amrex::Real m_fixed_roll{0};
    amrex::Real m_fixed_tilt{0};
    amrex::Real m_hub_yaw{0};
    amrex::Real m_hub_roll{0};
    amrex::Real m_hub_tilt{0};
    std::string m_spinner_mode{"fixed"};
    std::string m_turbine_label;
    bool m_hub_debug{false};
    vs::Vector m_hub_location{0, 0, 0};
    vs::Vector m_hub_location_init{0, 0, 0};

#ifdef AMR_WIND_USE_OPENFAST

    float turbine_base_pos[3]{0, 0, 0};
    float current_hub_abs_pos[3]{0, 0, 0};
    float current_hub_rot_vel[3]{0, 0, 0};
    double current_hub_orient[9]{0, 0, 0, 0, 0, 0, 0, 0, 0};

#endif

    amrex::Real m_last_hub_yaw{0};
    amrex::Real m_last_hub_roll{0};
    amrex::Real m_last_hub_tilt{0};
    amrex::Vector<amrex::Real> m_last_lidar_center{0.0, 0.0, 0.0};
    amrex::Vector<amrex::Real> step_lidar_center{0.0, 0.0, 0.0};
};

} // namespace amr_wind::sampling

#endif /* DTUSPINNERSAMPLER_H */
