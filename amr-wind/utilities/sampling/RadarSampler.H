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

    /** Read user inputs and initialize the sampling object
     *
     *  \param key Prefix used to parse inputs from file
     */
    void initialize(const std::string& key) override;

    enum class phase { FORWARD, FORWARD_PAUSE, REVERSE, REVERSE_PAUSE };
    //enum class NormalRule { SIGMA1, SIGMA2, SIGMA3, HALFPOWER };

    void check_bounds();

    void origin_cone();
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
    void output_locations(SampleLocType& locs) const override {return cone_axis_locations(locs);}

    void
    define_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;
    void
    populate_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;
    void output_netcdf_data(
        const ncutils::NCGroup& /*unused*/,
        const size_t /*unused*/) const override;

    //! Name of this sampling object
    std::string label() const override { return m_label; }
    std::string& label() override { return m_label; }

    //! Type of this sampling object
    std::string sampletype() const override { return identifier(); }

    //! Unique identifier for this set of probe locations
    int id() const override { return m_id; }
    int& id() override { return m_id; }

    //! Number of sampling container locations in a timestep
    int num_points() const override { return m_ntotal*(m_npts + ntheta*m_npts*(nphi-1)); }

    //! Number of probe locations to output
    int num_output_points() const override { return m_ntotal*m_npts; }

    //! Number of probe locations on Radar axis
    int num_points_axis() const { return m_npts; }

    //! Number of probe locations in a cone 
    int num_points_cone() const { return m_npts + ntheta*m_npts*(nphi-1); }

    //! Number of probe locations in a spherical cap    
    int num_points_quad() const { return 1 + ntheta*(nphi-1); }

    void
    line_average(
    const std::vector<double>& weights,
    const std::vector<double>& values,
    std::vector<double>& reduced,
    int offset);

    //! Run data modification for specific sampler
    virtual bool do_data_modification() { return true; }

    //! Modify sample buffer after sampling happens
    std::vector<double> modify_sample_data(const std::vector<double>&) override;




protected:
    const CFDSim& m_sim;

    std::string m_label;
    int m_id{-1};

    amrex::Real radtodeg{180.0 / M_PI};

    amrex::Vector<amrex::Real> m_start;
    amrex::Vector<amrex::Real> m_end;

    SamplerBase::SampleLocType initial_cone;
    SamplerBase::SampleLocType current_cones;

    std::vector<double> m_weights;
    std::vector<vs::Vector> m_rays;

    amrex::Real m_sample_freq; // Simulation sample rate
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

    amrex::Real m_radar_time; // Current radar time
    amrex::Real m_periodic_time; // Running cycle time (Forward, pause, reverse, pause)
    int m_ns;
    int m_ntotal;
    int m_cone_size;
    phase m_current_phase;

};

} // namespace amr_wind::sampling

#endif /* RADARSAMPLER_H */
