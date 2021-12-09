#include "amr-wind/utilities/sampling/DTUSpinnerSampler.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace sampling {

DTUSpinnerSampler::DTUSpinnerSampler(const CFDSim& sim) : LidarSampler(sim) {}

void DTUSpinnerSampler::initialize(const std::string& key)
{
    // Initialize the sampling time to be the same as simulation time
    //~ m_time_sampling = m_sim.time().current_time();

    LidarSampler::initialize(key);

    //~ subsampling();


    //~ // Read in new inputs specific to this class
    //~ amrex::ParmParse pp(key);

    //~ // This is the origin of the scan (x, y, z) [m]
    //~ pp.getarr("origin", m_origin);
    //~ AMREX_ALWAYS_ASSERT(static_cast<int>(m_origin.size()) == AMREX_SPACEDIM);

    //~ // The number of points
    //~ pp.get("num_points", m_npts);

    //~ // The length of the beam [m]
    //~ pp.get("length", m_length);

    //~ // The time table [s]
    //~ pp.getarr("time_table", m_time_table);
    //~ // Azimuth angle (this is in the N, E, S, W direction) [degrees]
    //~ pp.getarr("azimuth_table", m_azimuth_table);
    //~ // Elevation angle [degrees]
    //~ pp.getarr("elevation_table", m_elevation_table);

    //~ // Number of elements in the table
    //~ int np = m_time_table.size();

    //~ // Ensure that the tables have the same size
    //~ if (m_azimuth_table.size() != np) {
        //~ amrex::Abort(
            //~ "azimuth_table must have same number of entries as time_table ");
    //~ }
    //~ if (m_elevation_table.size() != np) {
        //~ amrex::Abort(
            //~ "elevation_table must have same number of entries as time_table ");
    //~ }

    //~ DTUSpinnerSampler::update_sampling_locations();

    //~ check_bounds();
}




void DTUSpinnerSampler::subsampling()
{
    amrex::Real time = m_sim.time().current_time();
    amrex::Real dt_sim = m_sim.time().deltaT();

    // Determine the number of subsamples
    m_ns = 0;

    // Loop to see how many times we will subsample
    while (m_time_sampling  < time + dt_sim) {
//~ std::cout << "m_time_sampling " << m_time_sampling << " time + dt_sim " <<  time + dt_sim << std::endl;
        m_ns +=1;
        m_time_sampling += m_dt_s;
    }

    // Resize these variables so they can store all the locations
    m_start.resize(AMREX_SPACEDIM * m_ns);
    m_end.resize(AMREX_SPACEDIM * m_ns);

std::cout <<" m_ns " << m_ns << std::endl;    

}

void DTUSpinnerSampler::sampling_locations(SampleLocType& locs) const
{

std::cout <<" m_npts m_ns "<< m_npts << " " << m_ns << std::endl;    

    // Resize to number of points in line times number of sampling times
    locs.resize(m_npts * m_ns);

    const amrex::Real ndiv = amrex::max(m_npts - 1, 1);
    amrex::Array<amrex::Real, AMREX_SPACEDIM> dx;

    // Loop per subsampling
    for (int k = 0; k < m_ns; ++k)
    {
        // Loop per spacial dimension
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            dx[d] = (m_end[d + k * m_npts] - m_start[d + k * m_npts]) / ndiv;
    
        for (int i = 0; i < m_npts; ++i) {
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                locs[i + k * m_ns][d] = m_start[d + k * m_npts] + i * dx[d];
        }
        
    }
}





void DTUSpinnerSampler::update_sampling_locations()
{

    amrex::Real time = m_sim.time().current_time();
    
    // Loop per subsampling
    for (int k = 0; k < m_ns; ++k)
    {
        time += m_dt_s * k;
        
        // The current azimuth angle
        const amrex::Real current_azimuth =
            ::amr_wind::interp::linear(m_time_table, m_azimuth_table, time);
    
        const amrex::Real current_elevation =
            ::amr_wind::interp::linear(m_time_table, m_elevation_table, time);
    
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            // Need to assign start point as the origin
            m_start[d + k * m_npts] = m_origin[d];
            // Initialize the end point
            m_end[d + k * m_npts] = m_origin[d];
        }
    
        // End point of the beam
        vs::Vector beam_vector = {m_length, 0., 0.};
    
        // The rotation matrix (takes in angles in degrees)
        vs::Tensor r1 = vs::yrot(current_elevation) & vs::zrot(current_azimuth);
    
        // Perform the vector rotation
        beam_vector = r1 & beam_vector;
    
        // Add the origin location to the beam vector
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            beam_vector[d] += m_origin[d + k * m_npts];
            m_end[d + k * m_npts] = beam_vector[d];
        }
    }
}

#ifdef AMR_WIND_USE_NETCDF
void DTUSpinnerSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
{
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("start", m_start);
    grp.put_attr("end", m_end);
    grp.put_attr("time_table", m_time_table);
    grp.put_attr("azimuth_table", m_azimuth_table);
    grp.put_attr("elevation_table", m_elevation_table);
    grp.def_var("points", NC_DOUBLE, {"num_time_steps", "num_points", "ndim"});
}

void DTUSpinnerSampler::populate_netcdf_metadata(const ncutils::NCGroup&) const {}
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
void DTUSpinnerSampler::define_netcdf_metadata(const ncutils::NCGroup&) const {}
void DTUSpinnerSampler::populate_netcdf_metadata(const ncutils::NCGroup&) const {}
void DTUSpinnerSampler::output_netcdf_data(
    const ncutils::NCGroup&, const size_t) const
{}
#endif

} // namespace sampling

template struct ::amr_wind::sampling::SamplerBase::Register<
    ::amr_wind::sampling::DTUSpinnerSampler>;

} // namespace amr_wind
