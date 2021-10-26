#include "amr-wind/utilities/sampling/IsoLineSampler.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tensor_ops.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace sampling {

IsoLineSampler::IsoLineSampler(const CFDSim& sim) : m_sim(sim) {}

IsoLineSampler::~IsoLineSampler() = default;

void IsoLineSampler::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);

    pp.get("num_points", m_npts);
    pp.getarr("start", m_start);
    pp.getarr("end", m_end);
    pp.getarr("orientation", m_oris);
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
        amrex::Print() << "WARNING: IsoLineSampler: Out of domain line was "
                          "truncated to match domain"
                       << std::endl;
    }
}

void IsoLineSampler::sampling_locations(SampleLocType& locs) const
{
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

void IsoLineSampler::sampling_orientations(SampleLocType& oris) const
{
    oris.resize(m_npts);

    for (int i = 0; i < m_npts; ++i) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            // orientation is constant for all points in line
            oris[i][d] = m_oris[d];
    }
}

#ifdef AMR_WIND_USE_NETCDF
void IsoLineSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
{
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("start", m_start);
    grp.put_attr("end", m_end);
}

void IsoLineSampler::populate_netcdf_metadata(const ncutils::NCGroup&) const {}
#else
void IsoLineSampler::define_netcdf_metadata(const ncutils::NCGroup&) const {}
void IsoLineSampler::populate_netcdf_metadata(const ncutils::NCGroup&) const {}
#endif

} // namespace sampling
} // namespace amr_wind
