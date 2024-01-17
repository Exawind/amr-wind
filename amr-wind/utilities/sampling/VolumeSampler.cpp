#include <limits>

#include "amr-wind/utilities/sampling/VolumeSampler.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::sampling {

VolumeSampler::VolumeSampler(const CFDSim& /*unused*/) {}

VolumeSampler::~VolumeSampler() = default;

void VolumeSampler::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);

    pp.getarr("axis", m_axis);
    pp.getarr("origin", m_origin);
    pp.getarr("num_points", m_npts_dir);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_axis.size()) == AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_origin.size()) == AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_npts_dir.size()) == AMREX_SPACEDIM);

    // Update total number of points
    const size_t tmp = m_npts_dir[0] * m_npts_dir[1] * m_npts_dir[2];
    if (tmp > static_cast<size_t>(std::numeric_limits<int>::max())) {
        amrex::Abort(
            "VolumeSampler: Plane definition " + key +
            " exceeds 32-bit integer limits");
    }
    m_npts = static_cast<int>(tmp);
}

void VolumeSampler::sampling_locations(SampleLocType& locs) const
{
    locs.resize(m_npts);

    amrex::Array<amrex::Real, AMREX_SPACEDIM> dx;

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        dx[d] = m_axis[d] / m_npts_dir[d];
    }

    int idx = 0;
    for (int k = 0; k < m_npts_dir[2]; ++k) {
        for (int j = 0; j < m_npts_dir[1]; ++j) {
            for (int i = 0; i < m_npts_dir[0]; ++i) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    locs[idx][0] = m_origin[0] + dx[0] * i;
                    locs[idx][1] = m_origin[1] + dx[1] * j;
                    locs[idx][2] = m_origin[2] + dx[2] * k;
                }
                ++idx;
            }
        }
    }
}

#ifdef AMR_WIND_USE_NETCDF
void VolumeSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
{
    const std::vector<int> ijk{m_npts_dir[0], m_npts_dir[1], m_npts_dir[2]};
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("ijk_dims", ijk);
    grp.put_attr("origin", m_origin);
    grp.put_attr("axis", m_axis);
}

void VolumeSampler::populate_netcdf_metadata(const ncutils::NCGroup&) const {}
#else
void VolumeSampler::define_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
void VolumeSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
#endif

} // namespace amr_wind::sampling
