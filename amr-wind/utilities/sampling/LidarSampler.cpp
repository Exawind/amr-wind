#include "amr-wind/utilities/sampling/LidarSampler.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace sampling {

void LidarSampler::initialize(const std::string& key)
{
	// Run the initialize function from the parent class
	LineSampler::initialize(key);

	// Read in new inputs specific to this class
    amrex::ParmParse pp(key);

	// This is the origin of the scan (x, y, z) [m]
	pp.getarr("origin", m_origin);

	// The length of the beam [m]
	pp.get("length", m_length);

	// This is the distance between points starting at zero [m]
	//~ pp.getarr("distribution", m_distribution);
	
	// The time table [s]
	pp.getarr("time_table", m_time_table);
	// Azimuth angle (this is in the N, E, S, W direction) [degrees]
	pp.getarr("azimuth_table", m_azimuth_table);
	// Elevation angle [degrees]
	pp.getarr("elevation_table", m_elevation_table);
	
	// Number of elements in the table
	m_npts = m_time_table.size();

	// Ensure that the tables have the same size
    if (m_azimuth_table.size() !=  m_npts) {
		amrex::Abort("azimuth_table must have same number of entries as time_table ");
    }
    if (m_elevation_table.size() !=  m_npts) {
		amrex::Abort("elevation_table must have same number of entries as time_table ");
    }

}


void LidarSampler::update_sampling_locations(SampleLocType& locs)
{
	amrex::Real time = m_sim.time().current_time();

	// The current azimuth angle
	m_current_azimuth = ::amr_wind::interp::linear(
            m_time_table, m_azimuth_table, time);

	m_current_elevation = ::amr_wind::interp::linear(
            m_time_table, m_elevation_table, time);


	// Need to assign start point and end point based on angles
	m_start = m_origin;

	//~ amrex::Real dx = m_length / m_npts;

	// End point of the beam 
	vs::Vector beam_vector = {m_length, 0., 0.};	
	
	// The rotation matrix (takes in angles in degrees)
	vs::Tensor r1 = vs::yrot(m_current_elevation) & vs::zrot(m_current_azimuth);
	
	// Perform the vector rotation
    beam_vector = r1 & beam_vector;

	// Add the origin location to the beam vector
	for (int d = 0; d < AMREX_SPACEDIM; ++d)
	{
        beam_vector[d] += m_origin[d];
        m_end[d] = beam_vector[d];
	}

	// Stup the locations from stat to end
	sampling_locations(locs);
	
}

} // namespace sampling
} // namespace amr_wind
