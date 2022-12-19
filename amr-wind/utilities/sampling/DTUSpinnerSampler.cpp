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
    m_time_sampling = m_sim.time().current_time();

    /*
     * Use this as a guide to implement the reading of the input
     * variables for the scanning pattern
     */

    // Read in new inputs specific to this class
    amrex::ParmParse pp(key);

    // This is the origin of the scan (x, y, z) [m]
    pp.getarr("origin", m_origin);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_origin.size()) == AMREX_SPACEDIM);

    // The number of points
    pp.get("num_points", m_npts);

    // The time step of the sampling
    pp.query("dt_s", m_dt_s);

    // The length of the beam [m]
    pp.get("length", m_length);

    // The time table [s]
    pp.getarr("time_table", m_time_table);
    // Azimuth angle (this is in the N, E, S, W direction) [degrees]
    pp.getarr("azimuth_table", m_azimuth_table);
    // Elevation angle [degrees]
    pp.getarr("elevation_table", m_elevation_table);

    // Number of elements in the table
    int np = static_cast<int>(m_time_table.size());

    // Ensure that the tables have the same size
    if (m_azimuth_table.size() != np) {
        amrex::Abort(
            "azimuth_table must have same number of entries as time_table ");
    }
    if (m_elevation_table.size() != np) {
        amrex::Abort(
            "elevation_table must have same number of entries as time_table ");
    }

    update_sampling_locations();

    check_bounds();
}

void DTUSpinnerSampler::sampling_locations(SampleLocType& locs) const
{

    // The total number of points at this time step
    int n_samples = m_npts * m_ns;

    // Resize to number of points in line times number of sampling times
    if (locs.size() < n_samples) {
        locs.resize(n_samples);
    }

    const amrex::Real ndiv = amrex::max(m_npts - 1, 1);
    amrex::Array<amrex::Real, AMREX_SPACEDIM> dx;

    // Loop per subsampling
    for (int k = 0; k < m_ns; ++k) {

        int offset = k * AMREX_SPACEDIM;

        // Loop per spacial dimension
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            dx[d] = (m_end[d + offset] - m_start[d + offset]) / ndiv;
        }

        for (int i = 0; i < m_npts; ++i) {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                locs[i + k * m_npts][d] = m_start[d + offset] + i * dx[d];
            }
        }
    }
}

void DTUSpinnerSampler::update_sampling_locations()
{

    BL_PROFILE(
        "amr-wind::Sampling::DTUSpinnerSampler::update_sampling_locations");

    amrex::Real time = m_sim.time().current_time();
    amrex::Real start_time = m_sim.time().start_time();
    amrex::Real dt_sim = m_sim.time().deltaT();
    // Determine the number of subsamples
    m_ns = 0;

    // Initialize the sampling time to the first time in the simulation
    if (time == start_time) {
        m_time_sampling = time;
    }
    amrex::Real time_tmp = m_time_sampling;
    bool cond = true;
    constexpr double eps = 1.0e-12;
    amrex::Real time_new = time + dt_sim;
    // Loop to see how many times we will subsample
    while (cond) {
        m_ns += 1;
        time_tmp += m_dt_s;
        cond = ((time_tmp + eps) < time_new);
    }

    int n_size = AMREX_SPACEDIM * m_ns;
    // Resize these variables so they can store all the locations
    if (m_start.size() < n_size) {
        m_start.resize(n_size);
    }
    if (m_end.size() < n_size) {
        m_end.resize(n_size);
    }

    // Loop per subsampling
    for (int k = 0; k < m_ns; ++k) {

        int offset = k * AMREX_SPACEDIM;

        m_time_sampling += m_dt_s;

        // The current azimuth angle
        const amrex::Real current_azimuth = ::amr_wind::interp::linear(
            m_time_table, m_azimuth_table, m_time_sampling);

        const amrex::Real current_elevation = ::amr_wind::interp::linear(
            m_time_table, m_elevation_table, m_time_sampling);

        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            // Need to assign start point as the origin
            m_start[d + offset] = m_origin[d];
            // Initialize the end point
            m_end[d + offset] = m_origin[d];
        }

        // End point of the beam
        vs::Vector beam_vector = {m_length, 0., 0.};

        // The rotation matrix (takes in angles in degrees)
        vs::Tensor r1 = vs::yrot(current_elevation) & vs::zrot(current_azimuth);

        // Perform the vector rotation
        beam_vector = r1 & beam_vector;

        // Add the origin location to the beam vector
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {

            beam_vector[d] += m_origin[d];
            m_end[d + offset] = beam_vector[d];
        }
    }
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
