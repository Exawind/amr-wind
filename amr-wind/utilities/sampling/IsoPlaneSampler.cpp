#include <limits>

#include "amr-wind/utilities/sampling/IsoPlaneSampler.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace sampling {

IsoPlaneSampler::IsoPlaneSampler(const CFDSim&) {}

IsoPlaneSampler::~IsoPlaneSampler() = default;

void IsoPlaneSampler::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);

    pp.getarr("axis1", m_axis1);
    pp.getarr("axis2", m_axis2);
    pp.getarr("origin", m_origin);
    pp.getarr("num_points", m_npts_dir);
    pp.getarr("orientation", m_oris);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_axis1.size()) == AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_axis2.size()) == AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_origin.size()) == AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_npts_dir.size()) == 2);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_oris.size()) == AMREX_SPACEDIM);

    // Normalize orientation
    amrex::Real mag = 0;
    for (int i = 0; i < m_oris.size(); ++i) {
        // Get norm
        mag += std::pow(m_oris[i], 2);
    }
    mag = std::sqrt(mag);
    for (int i = 0; i < m_oris.size(); ++i) {
        // Normalize orientation vector
        m_oris[i] /= mag;
    }

    // Update total number of points
    const size_t tmp = m_npts_dir[0] * m_npts_dir[1];
    if (tmp > static_cast<size_t>(std::numeric_limits<int>::max())) {
        amrex::Abort(
            "IsoPlaneSampler: Plane definition " + key +
            " exceeds 32-bit integer limits");
    }
    m_npts = static_cast<int>(tmp);
}

void IsoPlaneSampler::sampling_locations(SampleLocType& locs) const
{
    locs.resize(m_npts);

    amrex::Array<amrex::Real, AMREX_SPACEDIM> dx;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> dy;

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        dx[d] = m_axis1[d] / amrex::max(m_npts_dir[0] - 1, 1);
        dy[d] = m_axis2[d] / amrex::max(m_npts_dir[1] - 1, 1);
    }

    int idx = 0;
    for (int j = 0; j < m_npts_dir[1]; ++j) {
        for (int i = 0; i < m_npts_dir[0]; ++i) {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                locs[idx][d] = m_origin[d] + dx[d] * i + dy[d] * j;
            }
            ++idx;
        }
    }
}

void IsoPlaneSampler::sampling_orientations(SampleLocType& oris) const
{
    oris.resize(m_npts);

    int idx = 0;
    for (int j = 0; j < m_npts_dir[1]; ++j) {
        for (int i = 0; i < m_npts_dir[0]; ++i) {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                oris[idx][d] = m_oris[d];
            }
            ++idx;
        }
    }
}

#ifdef AMR_WIND_USE_NETCDF
void IsoPlaneSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
{
    const std::vector<int> ijk{
        m_npts_dir[0], m_npts_dir[1], static_cast<int>(m_poffsets.size())};
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("ijk_dims", ijk);
    grp.put_attr("origin", m_origin);
    grp.put_attr("axis1", m_axis1);
    grp.put_attr("axis2", m_axis2);
    grp.put_attr("orientation", m_oris);
}

void IsoPlaneSampler::populate_netcdf_metadata(const ncutils::NCGroup&) const {}
#else
void IsoPlaneSampler::define_netcdf_metadata(const ncutils::NCGroup&) const {}
void IsoPlaneSampler::populate_netcdf_metadata(const ncutils::NCGroup&) const {}
#endif

} // namespace sampling
} // namespace amr_wind
