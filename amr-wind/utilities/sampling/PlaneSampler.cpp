#include "amr-wind/utilities/sampling/PlaneSampler.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace sampling {

PlaneSampler::PlaneSampler(const CFDSim&) {}

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
    int noffsets = m_poffsets.size();
    if (noffsets > 0) {
        pp.getarr("normal", m_normal);
        AMREX_ALWAYS_ASSERT(
            static_cast<int>(m_normal.size()) == AMREX_SPACEDIM);
    } else {
        m_poffsets.push_back(0.0);
    }

    // Update total number of points
    m_npts = m_npts_dir[0] * m_npts_dir[1] * m_poffsets.size();
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
    const int nplanes = m_poffsets.size();
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

} // namespace sampling
} // namespace amr_wind
