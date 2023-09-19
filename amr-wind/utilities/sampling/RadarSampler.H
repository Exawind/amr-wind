#ifndef RADARSAMPLER_H
#define RADARSAMPLER_H

#include "amr-wind/utilities/sampling/SamplerBase.H"
#include "amr-wind/utilities/sampling/SamplingUtils.H"
#include "amr-wind/utilities/sampling/LineSampler.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind::sampling {

/** Sample data along a line
 *  \ingroup sampling
 *
 *  Defines probes along a line defined by `start` and `end` coordinates and
 *  divides it into equal segements as determined by `num_points`. The length of
 *  each segment is given by \f[ (end - start) / (num_points - 1) \f]
 */
class RadarSampler : public SamplerBase::Register<RadarSampler>
{
public:
    static std::string identifier() { return "RadarSampler"; }

    explicit RadarSampler(const CFDSim& /*sim*/);

    ~RadarSampler() override;

    using LosRotType = amrex::Vector<vs::Tensor>;
    using LosUnitType = amrex::Vector<vs::Vector>;

    /** Read user inputs and initialize the sampling object
     *
     *  \param key Prefix used to parse inputs from file
     */
    void initialize(const std::string& key) override;

    enum class phase { FORWARD, FORWARD_PAUSE, REVERSE, REVERSE_PAUSE };

    void check_bounds();

    void new_cone();

    double total_sweep_time() const;
    double periodic_time();
    int sweep_count() const;
    phase determine_operation_phase() const;
    double determine_current_sweep_angle() const;

    //! Populate and return a vector of probe locations to be sampled
    void sampling_locations(SampleLocType& /*locs*/) const override;
    void update_sampling_locations() override;
    void cone_axis_locations(SampleLocType& /*locs*/) const;
    void output_locations(SampleLocType& locs) const override
    {
        return cone_axis_locations(locs);
    }

    void post_sample_actions() override;

    void
    define_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;
    void
    populate_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;
    void output_netcdf_data(
        const ncutils::NCGroup& /*unused*/,
        const size_t /*unused*/) const override;

    bool output_netcdf_field(
        const std::vector<double>& output_buffer,
        ncutils::NCGroup& grp,
        const size_t nt) override;

    //! Name of this sampling object
    std::string label() const override { return m_label; }
    std::string& label() override { return m_label; }

    //! Type of this sampling object
    std::string sampletype() const override { return identifier(); }

    //! Unique identifier for this set of probe locations
    int id() const override { return m_id; }
    int& id() override { return m_id; }

    //! Number of sampling container locations in a timestep
    long num_points() const override
    {
        return m_interp_steps * m_ntotal *
               (m_npts + ntheta * m_npts * (nphi - 1));
    }

    //! Number of probe locations to output
    long num_output_points() const override { return m_ntotal * m_npts; }

    //! Number of sampling container locations per interp step
    long num_points_scan() const
    {
        return m_ntotal * (m_npts + ntheta * m_npts * (nphi - 1));
    }

    //! Number of probe locations on Radar axis
    long num_points_axis() const { return m_npts; }

    //! Number of probe locations in a cone
    long num_points_cone() const
    {
        return m_npts + ntheta * m_npts * (nphi - 1);
    }

    //! Number of probe locations in a spherical cap
    long num_points_quad() const { return 1 + ntheta * (nphi - 1); }

    void line_average(
        const std::vector<double>& weights,
        const std::vector<double>& values,
        std::vector<double>& reduced,
        int offset);

    //! Run data modification for specific sampler
    bool do_data_modification() override { return true; }
    bool do_convert_velocity_los() override { return true; }
    bool do_subsampling_interp() override { return true; }

    //! Modify sample buffer after sampling happens
    std::vector<double> modify_sample_data(
        const std::vector<double>& sample_data,
        const std::string& /*unused*/) override;

    void calc_lineofsight_velocity(
        const std::vector<std::vector<double>>&, const int) override;

protected:
    const CFDSim& m_sim;

    std::string m_label;
    int m_id{-1};

    amrex::Real radtodeg{180.0 / M_PI};

    amrex::Vector<amrex::Real> m_start;
    amrex::Vector<amrex::Real> m_end;

    SamplerBase::SampleLocType initial_cone;
    SamplerBase::SampleLocType current_cones;
    SamplerBase::SampleLocType prior_cones;
    SamplerBase::SampleLocType sample_cones;

    LosRotType los_proj;
    LosUnitType los_unit;

    std::vector<double> m_los_velocity;
    std::vector<double> m_los_velocity_next;
    std::vector<double> m_los_velocity_prior;
    std::vector<double> m_los_velocity_interp;

    std::vector<double> m_weights;
    std::vector<vs::Vector> m_rays;

    amrex::Real m_sample_freq;       // Simulation sample rate
    amrex::Real m_radar_sample_freq; // Actual device sample rate
    int m_npts{0};
    amrex::Real m_cone_angle;
    std::string m_radar_quad_type;
    int m_npts_azimuth{0};
    amrex::Real m_beam_length;
    amrex::Real m_angular_speed;
    amrex::Vector<amrex::Real> m_axis;
    amrex::Vector<amrex::Real> m_vertical;
    amrex::Real m_sweep_angle;
    amrex::Real m_reset_time;
    amrex::Vector<amrex::Real> m_elevation_angles;

    int nphi{5};
    int ntheta{5};

    vs::Vector m_radar_rotation_axis;

    long m_radar_iter{0};     // Iteration tracking from radar perspective
    amrex::Real m_radar_time; // Current radar time
    amrex::Real
        m_periodic_time; // Running cycle time (Forward, pause, reverse, pause)
    int m_ns;            // Number of radar steps in current timestep
    int m_ntotal;        // Max number of radar timesteps in any timestep
    int m_cone_size;
    phase m_current_phase;

    bool m_debug_print{false};
    bool m_output_cone_points{false};

    long m_interp_steps{2};
    int m_ns_prior;
    int m_ntotal_prior;
    long m_num_points_prior;
    long m_num_output_points_prior;

    double m_fill_val{-99999.99};
};

} // namespace amr_wind::sampling

#endif /* RADARSAMPLER_H */
