#ifndef LIDARSAMPLER_H
#define LIDARSAMPLER_H

//~ #include "amr-wind/utilities/sampling/SamplerBase.H"
#include "amr-wind/utilities/sampling/LineSampler.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind {
namespace sampling {

/** Sample data along a line
 *  \ingroup sampling
 *
 *  Defines probes along a line defined by `start` and `end` coordinates and
 *  divides it into equal segements as determined by `num_points`. The length of
 *  each segment is given by \f[ (end - start) / (num_points - 1) \f]
 */
//~ class LidarSampler : public SamplerBase::Register<LidarSampler>
//~ class LidarSampler : public LineSampler::Register<LidarSampler>
class LidarSampler : public LineSampler
{
public:
    static const std::string identifier() { return "LidarSampler"; }

    //~ explicit LidarSampler(const CFDSim&);

    //~ virtual ~LidarSampler();

    //~ /** Read user inputs and initialize the sampling object
     //~ *
     //~ *  \param key Prefix used to parse inputs from file
     //~ */
    void initialize(const std::string& key) override;


	/** Update the position of the sampling locations in time
	 * 
	 * 
	 */
    void update_sampling_locations(SampleLocType& locs) override;

    //~ void update_positions();

    //~ //! Populate and return a vector of probe locations to be sampled
    //~ void sampling_locations(SampleLocType&) const override;

    //~ void define_netcdf_metadata(const ncutils::NCGroup&) const override;
    //~ void populate_netcdf_metadata(const ncutils::NCGroup&) const override;

    //~ //! Name of this sampling object
    //~ std::string label() const override { return m_label; }
    //~ std::string& label() override { return m_label; }

    //~ //! Unique identifier for this set of probe locations
    //~ int id() const override { return m_id; }
    //~ int& id() override { return m_id; }

    //~ //! Number of probe locations along the line
    //~ int num_points() const override { return m_npts; }

//~ private:
    //~ const CFDSim& m_sim;

    amrex::Vector<amrex::Real> m_origin;
    amrex::Vector<amrex::Real> m_time_table;
    amrex::Vector<amrex::Real> m_azimuth_table;
    amrex::Vector<amrex::Real> m_elevation_table;
    amrex::Vector<amrex::Real> m_distribution;

    //~ amrex::Vector<amrex::Real> m_start;
    //~ amrex::Vector<amrex::Real> m_end;

	amrex::Real m_current_azimuth{0};
	amrex::Real m_current_elevation{0};
	amrex::Real m_length{0};

    //~ std::string m_label;
    //~ int m_id{-1};
    //~ int m_npts{0};
};

} // namespace sampling
} // namespace amr_wind

#endif /* LINESAMPLER_H */
