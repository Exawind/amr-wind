#include "amr-wind/utilities/sampling/ConeSampler.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tensor_ops.H"

#include "AMReX_ParmParse.H"

#define PI 3.14159265

namespace amr_wind::sampling {

ConeSampler::ConeSampler(const CFDSim& sim) : m_sim(sim) {}

ConeSampler::~ConeSampler() = default;

void ConeSampler::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);

    pp.get("num_points", m_npts);
    pp.get("num_theta", m_ntheta);
    pp.get("num_phi", m_nphi);
    pp.get("phi", m_phi);
    pp.getarr("start", m_start);
    pp.getarr("end", m_end);

    check_bounds();
}

vs::Vector canon_rotation(const vs::Vector& angles, const vs::Vector& data)
{
    // Note: Angles in degrees as per "vs" class
    const vs::Tensor rotMatrix =
        vs::xrot(angles.x()) & vs::yrot(angles.y()) & vs::zrot(angles.z());
    return data & rotMatrix;
}

auto rotate_euler_vec(vs::Vector axis, double angle, vs::Vector vec)
{
    axis.normalize();
    const auto RotMat = vs::quaternion(axis, angle);
    return vec & RotMat;
}

void ConeSampler::check_bounds()
{
    const int lev = 0;
    const auto* prob_lo = m_sim.mesh().Geom(lev).ProbLo();
    const auto* prob_hi = m_sim.mesh().Geom(lev).ProbHi();

    bool all_ok = true;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (m_start[d] < prob_lo[d]) {
            all_ok = false;
            m_start[d] = prob_lo[d];
        }
        if (m_start[d] > prob_hi[d]) {
            all_ok = false;
            m_start[d] = prob_lo[d];
        }
        if (m_end[d] < prob_lo[d]) {
            all_ok = false;
            m_end[d] = prob_lo[d];
        }
        if (m_end[d] > prob_hi[d]) {
            all_ok = false;
            m_end[d] = prob_lo[d];
        }
    }
    if (!all_ok) {
        amrex::Print() << "WARNING: ConeSampler: Out of domain line was "
                          "truncated to match domain"
                       << std::endl;
    }
}

void ConeSampler::sampling_locations(SampleLocType& locs) const
{
    // size assumes only one axis set of points
    int full_size = m_npts + m_ntheta * (m_npts - 1) * (m_nphi - 1);
    locs.resize(full_size);

    vs::Vector canon_vector(0.0, 0.0, 1.0);
    vs::Vector origin_vector(0.0, 0.0, 0.0);

    const amrex::Real nptsdiv = amrex::max(m_npts - 1, 1);
    const amrex::Real nthediv = amrex::max(m_ntheta, 1);
    const amrex::Real nphidiv = amrex::max(m_nphi - 1, 1);
    amrex::Array<amrex::Real, AMREX_SPACEDIM> dx;

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        dx[d] = (canon_vector[d] - origin_vector[d]) / nptsdiv;
    }

    // Assemble cone axis unit beam
    for (int i = 0; i < m_npts; ++i) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            locs[i][d] = origin_vector[d] + i * dx[d];
        }
    }

    const amrex::Real dphi = m_phi / nphidiv;
    const amrex::Real dtheta = 360.0 / nthediv;
    amrex::Real c_theta = 0.0;
    amrex::Real c_phi = 0.0;

    // Assemble all other unit beams
    for (int t = 0; t < m_ntheta; ++t) {
        c_theta = t * dtheta;
        for (int p = 0; p < m_nphi - 1; ++p) {
            c_phi = static_cast<double>((p + 1) * dphi);
            for (int i = 0; i < m_npts - 1; ++i) {
                int ci = m_npts + i + (m_npts - 1) * (p + (m_nphi - 1) * t);

                const vs::Vector temp_loc(
                    origin_vector[0] + dx[0] + i * dx[0],
                    origin_vector[1] + dx[1] + i * dx[1],
                    origin_vector[2] + dx[2] + i * dx[2]);
                const vs::Vector temp_ang(0.0, c_phi, c_theta);
                vs::Vector rot_loc = canon_rotation(temp_ang, temp_loc);

                locs[ci][0] = rot_loc[0];
                locs[ci][1] = rot_loc[1];
                locs[ci][2] = rot_loc[2];
            }
        }
    }

    // Transform to scale and align cone axis with start and end points
    vs::Vector target_beam(
        m_end[0] - m_start[0], m_end[1] - m_start[1], m_end[2] - m_start[2]);
    amrex::Real target_beam_scale = std::sqrt(
        std::pow(target_beam[0], 2.0) + std::pow(target_beam[1], 2.0) +
        std::pow(target_beam[2], 2.0));

    vs::Vector unit_target = target_beam.normalize();
    vs::Vector tc_axis = canon_vector ^ unit_target;
    amrex::Real tc_angle = std::acos(unit_target & canon_vector) * 180.0 / PI;

    for (int i = 0; i < full_size; ++i) {
        vs::Vector temp_loc(locs[i][0], locs[i][1], locs[i][2]);
        vs::Vector new_rot(rotate_euler_vec(tc_axis, tc_angle, temp_loc));
        locs[i][0] = new_rot[0] * target_beam_scale + m_start[0];
        locs[i][1] = new_rot[1] * target_beam_scale + m_start[1];
        locs[i][2] = new_rot[2] * target_beam_scale + m_start[2];
    }
}

#ifdef AMR_WIND_USE_NETCDF
void ConeSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
{
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("start", m_start);
    grp.put_attr("end", m_end);

    grp.def_dim("ndim", AMREX_SPACEDIM);
    grp.def_var("points", NC_DOUBLE, {"num_time_steps", "num_points", "ndim"});
}

void ConeSampler::populate_netcdf_metadata(const ncutils::NCGroup&) const {}

void ConeSampler::output_netcdf_data(
    const ncutils::NCGroup& grp, const size_t nt) const
{
    amrex::Print() << "Output Cone Sampler Netcdf Data" << std::endl;
    std::vector<size_t> start{nt, 0, 0};
    std::vector<size_t> count{1, 0, AMREX_SPACEDIM};
    SamplerBase::SampleLocType locs;
    sampling_locations(locs);

    auto xyz = grp.var("points");
    count[1] = num_points();
    xyz.put(&locs[0][0], start, count);
}

#else

void ConeSampler::define_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
void ConeSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
void ConeSampler::output_netcdf_data(
    const ncutils::NCGroup& /*unused*/, const size_t /*unused*/) const
{}

#endif

} // namespace amr_wind::sampling
