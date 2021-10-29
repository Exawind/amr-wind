#include "amr-wind/utilities/sampling/LidarSampler.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace sampling {

//~ LidarSampler::LidarSampler(const CFDSim& sim) : m_sim(sim) {}

//~ LidarSampler::~LidarSampler() = default;

void LidarSampler::initialize(const std::string& key)
{
	// Run the initialize function from the parent class
	LineSampler::initialize(key);

	// Read in new inputs specific to this class
    amrex::ParmParse pp(key);

	// This vector points the direction of the scan
	// Will be converted to non-dimensional
	//~ pp.getarr("orientation_vector", m_ov);
	//** need to implement azimuth and elevation as function of time [angles in degrees]
	

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
	// Elevation angle
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


    //~ if (m_npts == 1) {
        //~ m_current_azimuth = meta.thrust_coeff[0];
        //~ meta.table_velocity = {1.0};
    //~ }

    //~ pp.get("num_points", m_npts);
    //~ pp.getarr("start", m_start);
    //~ pp.getarr("end", m_end);

    //~ const int lev = 0;
    //~ const auto* prob_lo = m_sim.mesh().Geom(lev).ProbLo();
    //~ const auto* prob_hi = m_sim.mesh().Geom(lev).ProbHi();

    //~ bool all_ok = true;
    //~ for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        //~ if (m_start[d] < prob_lo[d]) {
            //~ all_ok = false;
            //~ m_start[d] = prob_lo[d];
        //~ }
        //~ if (m_start[d] > prob_hi[d]) {
            //~ all_ok = false;
            //~ m_start[d] = prob_lo[d];
        //~ }
        //~ if (m_end[d] < prob_lo[d]) {
            //~ all_ok = false;
            //~ m_end[d] = prob_lo[d];
        //~ }
        //~ if (m_end[d] > prob_hi[d]) {
            //~ all_ok = false;
            //~ m_end[d] = prob_lo[d];
        //~ }
    //~ }
    //~ if (!all_ok) {
        //~ amrex::Print() << "WARNING: LineSampler: Out of domain line was "
                          //~ "truncated to match domain"
                       //~ << std::endl;
    //~ }
}



//~ void LidarSampler::update_positions()
//~ {
//~ }



void LidarSampler::update_sampling_locations(SampleLocType& locs)
{
	// Update the positions of all 
    //~ amrex::ParmParse pp(key);

	//~ const CFDSim& m_sim;
	//~ const auto& time = LineSampler::m_sim.time();
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

    locs.resize(m_npts);

    const amrex::Real ndiv = amrex::max(m_npts - 1, 1);
    amrex::Array<amrex::Real, AMREX_SPACEDIM> dx;

    for (int d = 0; d < AMREX_SPACEDIM; ++d)
        dx[d] = (m_end[d] - m_start[d]) / ndiv;

    for (int i = 0; i < m_npts; ++i) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            locs[i][d] = m_start[d] + i * dx[d];
    }
	
}

//~ void LidarSampler::sampling_locations(SampleLocType& locs) const
//~ {
    //~ locs.resize(m_npts);

    //~ const amrex::Real ndiv = amrex::max(m_npts - 1, 1);
    //~ amrex::Array<amrex::Real, AMREX_SPACEDIM> dx;

    //~ for (int d = 0; d < AMREX_SPACEDIM; ++d)
        //~ dx[d] = (m_end[d] - m_start[d]) / ndiv;

    //~ for (int i = 0; i < m_npts; ++i) {
        //~ for (int d = 0; d < AMREX_SPACEDIM; ++d)
            //~ locs[i][d] = m_start[d] + i * dx[d];
    //~ }
//~ }

//~ #ifdef AMR_WIND_USE_NETCDF
//~ void LidarSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
//~ {
    //~ grp.put_attr("sampling_type", identifier());
    //~ grp.put_attr("start", m_start);
    //~ grp.put_attr("end", m_end);
//~ }

//~ void LidarSampler::populate_netcdf_metadata(const ncutils::NCGroup&) const {}
//~ #else
//~ void LidarSampler::define_netcdf_metadata(const ncutils::NCGroup&) const {}
//~ void LidarSampler::populate_netcdf_metadata(const ncutils::NCGroup&) const {}
//~ #endif

} // namespace sampling
} // namespace amr_wind
