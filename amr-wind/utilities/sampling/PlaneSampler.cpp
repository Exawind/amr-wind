#include <limits>

#include "amr-wind/utilities/sampling/PlaneSampler.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::sampling {

PlaneSampler::PlaneSampler(const CFDSim& /*unused*/) {}

PlaneSampler::~PlaneSampler() = default;

void PlaneSampler::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);

    pp.getarr("axis1", m_axis1);
    pp.getarr("axis2", m_axis2);
    pp.getarr("origin", m_origin);
    pp.getarr("num_points", m_npts_dir);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_axis1.size()) == AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_axis2.size()) == AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_origin.size()) == AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_npts_dir.size()) == 2);

    pp.queryarr("offsets", m_poffsets);
    int noffsets = static_cast<int>(m_poffsets.size());
    if (noffsets > 0) {
        pp.getarr("normal", m_normal);
        AMREX_ALWAYS_ASSERT(
            static_cast<int>(m_normal.size()) == AMREX_SPACEDIM);
    } else {
        m_poffsets.push_back(0.0);
    }

    // Update total number of points
    const size_t tmp = m_poffsets.size() * m_npts_dir[0] * m_npts_dir[1];
    if (tmp > static_cast<size_t>(std::numeric_limits<int>::max())) {
        amrex::Abort(
            "PlaneSampler: Plane definition " + key +
            " exceeds 32-bit integer limits");
    }
    m_npts = static_cast<int>(tmp);
}

void PlaneSampler::sampling_locations(SampleLocType& locs) const
{
    locs.resize(m_npts);

    amrex::Array<amrex::Real, AMREX_SPACEDIM> dx;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> dy;

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        dx[d] = m_axis1[d] / amrex::max(m_npts_dir[0] - 1, 1);
        dy[d] = m_axis2[d] / amrex::max(m_npts_dir[1] - 1, 1);
    }

    int idx = 0;
    const int nplanes = static_cast<int>(m_poffsets.size());
    for (int k = 0; k < nplanes; ++k) {
        for (int j = 0; j < m_npts_dir[1]; ++j) {
            for (int i = 0; i < m_npts_dir[0]; ++i) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    locs[idx][d] = m_origin[d] + dx[d] * i + dy[d] * j +
                                   m_poffsets[k] * m_normal[d];
                }
                ++idx;
            }
        }
    }
}

#ifdef AMR_WIND_USE_NETCDF
void PlaneSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
{
    const std::vector<int> ijk{
        m_npts_dir[0], m_npts_dir[1], static_cast<int>(m_poffsets.size())};
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("ijk_dims", ijk);
    grp.put_attr("origin", m_origin);
    grp.put_attr("axis1", m_axis1);
    grp.put_attr("axis2", m_axis2);
    grp.put_attr("axis3", m_normal);
    grp.put_attr("offsets", m_poffsets);
}

void PlaneSampler::populate_netcdf_metadata(const ncutils::NCGroup&) const {}
#else
void PlaneSampler::define_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
void PlaneSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
#endif

} // namespace amr_wind::sampling
