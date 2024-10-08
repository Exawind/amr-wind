#include "amr-wind/utilities/sampling/ProbeSampler.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::sampling {

ProbeSampler::ProbeSampler(const CFDSim& sim) : m_sim(sim) {}

ProbeSampler::~ProbeSampler() = default;

void ProbeSampler::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);
    std::string pfile("probe_locations.txt");
    pp.query("probe_location_file", pfile);

    std::ifstream ifh(pfile, std::ios::in);
    if (!ifh.good()) {
        amrex::Abort("Cannot find probe location file: " + pfile);
    }

    pp.queryarr("offsets", m_poffsets);
    if (m_poffsets.size() > 0) {
        pp.getarr("offset_vector", m_offset_vector);
        AMREX_ALWAYS_ASSERT(
            static_cast<int>(m_offset_vector.size()) == AMREX_SPACEDIM);
    } else {
        // No offsets is implemented as 1 offset of 0.
        m_poffsets.push_back(0.0);
    }

    int npts_file = 0;
    ifh >> npts_file;
    ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    SampleLocType probes_file;
    probes_file.resize(npts_file);
    m_npts = m_poffsets.size() * npts_file;
    m_probes.resize(m_npts);
    // Read through points in file
    for (int i = 0; i < npts_file; ++i) {
        ifh >> probes_file[i][0] >> probes_file[i][1] >> probes_file[i][2];
        ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    // Incorporate offsets
    for (int n = 0; n < m_poffsets.size(); ++n) {
        for (int i = 0; i < npts_file; ++i) {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                m_probes[i + n * npts_file][d] =
                    probes_file[i][d] + m_poffsets[n] * m_offset_vector[d];
            }
        }
    }

    check_bounds();
}

void ProbeSampler::check_bounds()
{
    const int lev = 0;
    const auto* prob_lo = m_sim.mesh().Geom(lev).ProbLo();
    const auto* prob_hi = m_sim.mesh().Geom(lev).ProbHi();
    bool all_ok = true;
    for (int i = 0; i < m_npts; ++i) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (m_probes[i][d] <= prob_lo[d]) {
                all_ok = false;
                m_probes[i][d] = prob_lo[d] + bounds_tol;
            }
            if (m_probes[i][d] >= prob_hi[d]) {
                all_ok = false;
                m_probes[i][d] = prob_hi[d] - bounds_tol;
            }
        }
    }
    if (!all_ok) {
        amrex::Print() << "WARNING: ProbeSampler: Out of domain probe was "
                          "truncated to match domain"
                       << std::endl;
    }
}

void ProbeSampler::sampling_locations(SampleLocType& locs) const
{
    locs.resize(m_npts);
    for (int i = 0; i < m_npts; ++i) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            locs[i][d] = m_probes[i][d];
        }
    }
}

#ifdef AMR_WIND_USE_NETCDF
void ProbeSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
{
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("offset_vector", m_offset_vector);
    grp.put_attr("offsets", m_poffsets);
}
#else
void ProbeSampler::define_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
#endif

} // namespace amr_wind::sampling
