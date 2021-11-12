#include "amr-wind/utilities/sampling/LidarSampler.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace sampling {

LidarSampler::LidarSampler(const CFDSim& sim) : LineSampler(sim) {}

void LidarSampler::initialize(const std::string& key)
{
    // Read in new inputs specific to this class
    amrex::ParmParse pp(key);

    // This is the origin of the scan (x, y, z) [m]
    pp.getarr("origin", m_origin);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_origin.size()) == AMREX_SPACEDIM);

    // The number of points
    pp.get("num_points", m_npts);

    // The length of the beam [m]
    pp.get("length", m_length);

    // The time table [s]
    pp.getarr("time_table", m_time_table);
    // Azimuth angle (this is in the N, E, S, W direction) [degrees]
    pp.getarr("azimuth_table", m_azimuth_table);
    // Elevation angle [degrees]
    pp.getarr("elevation_table", m_elevation_table);

    // Number of elements in the table
    int np = m_time_table.size();

    // Ensure that the tables have the same size
    if (m_azimuth_table.size() != np) {
        amrex::Abort(
            "azimuth_table must have same number of entries as time_table ");
    }
    if (m_elevation_table.size() != np) {
        amrex::Abort(
            "elevation_table must have same number of entries as time_table ");
    }

    LidarSampler::update_sampling_locations();

    check_bounds();
}

void LidarSampler::update_sampling_locations()
{

    amrex::Real time = m_sim.time().current_time();

    // The current azimuth angle
    const amrex::Real current_azimuth =
        ::amr_wind::interp::linear(m_time_table, m_azimuth_table, time);

    const amrex::Real current_elevation =
        ::amr_wind::interp::linear(m_time_table, m_elevation_table, time);

    // Need to assign start point as the origin
    m_start = m_origin;
    // Initialize the end point
    m_end = m_origin;

    // End point of the beam
    vs::Vector beam_vector = {m_length, 0., 0.};

    // The rotation matrix (takes in angles in degrees)
    vs::Tensor r1 = vs::yrot(current_elevation) & vs::zrot(current_azimuth);

    // Perform the vector rotation
    beam_vector = r1 & beam_vector;

    // Add the origin location to the beam vector
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        beam_vector[d] += m_origin[d];
        m_end[d] = beam_vector[d];
    }
}

#ifdef AMR_WIND_USE_NETCDF
void LidarSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
{
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("start", m_start);
    grp.put_attr("end", m_end);
    grp.put_attr("time_table", m_time_table);
    grp.put_attr("azimuth_table", m_azimuth_table);
    grp.put_attr("elevation_table", m_elevation_table);
    grp.def_var("points", NC_DOUBLE, {"num_time_steps", "num_points", "ndim"});
}

void LidarSampler::populate_netcdf_metadata(const ncutils::NCGroup&) const {}
void LidarSampler::output_netcdf_data(
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
void LidarSampler::define_netcdf_metadata(const ncutils::NCGroup&) const {}
void LidarSampler::populate_netcdf_metadata(const ncutils::NCGroup&) const {}
void LidarSampler::output_netcdf_data(
    const ncutils::NCGroup&, const size_t) const
{}
#endif

} // namespace sampling

template struct ::amr_wind::sampling::SamplerBase::Register<
    ::amr_wind::sampling::LidarSampler>;

} // namespace amr_wind
