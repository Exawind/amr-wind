#include "amr-wind/utilities/sampling/ProbeSampler.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace sampling {

ProbeSampler::ProbeSampler(const CFDSim&) {}

ProbeSampler::~ProbeSampler() = default;

void ProbeSampler::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);
    std::string pfile("probe_locations.txt");
    pp.query("probe_location_file", pfile);

    std::ifstream ifh(pfile, std::ios::in);
    if (!ifh.good()) amrex::Abort("Cannot find probe location file: " + pfile);

    ifh >> m_npts;
    ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    m_probes.resize(m_npts);
    for (int i = 0; i < m_npts; ++i) {
        ifh >> m_probes[i][0] >> m_probes[i][1] >> m_probes[i][2];
        ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

void ProbeSampler::initialize_iso(const std::string& key)
{
    amrex::ParmParse pp(key);
    std::string pfile("probe_locations.txt");
    pp.query("probe_location_file", pfile);
    std::string ofile("probe_orientations.txt");
    pp.query("probe_orientation_file", ofile);

    // Initial positions
    std::ifstream ifh(pfile, std::ios::in);
    if (!ifh.good()) amrex::Abort("Cannot find probe location file: " + pfile);

    ifh >> m_npts;
    ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    m_probes.resize(m_npts);
    for (int i = 0; i < m_npts; ++i) {
        ifh >> m_probes[i][0] >> m_probes[i][1] >> m_probes[i][2];
        ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    // Search directions (orientation)
    std::ifstream ifo(ofile, std::ios::in);
    if (!ifo.good())
        amrex::Abort("Cannot find probe orientation file: " + ofile);

    int m_opts;
    ifo >> m_opts;
    if (m_opts != m_npts)
        amrex::Abort("Probe location number conflicts with orientations");

    ifo.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    m_oris.resize(m_npts);
    for (int i = 0; i < m_npts; ++i) {
        ifo >> m_oris[i][0] >> m_oris[i][1] >> m_oris[i][2];
        ifo.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        // Normalize values
        amrex::Real mag = std::sqrt(
            pow(m_oris[i][0], 2) + pow(m_oris[i][1], 2) + pow(m_oris[i][2], 2));
        m_oris[i][0] /= mag;
        m_oris[i][1] /= mag;
        m_oris[i][2] /= mag;
    }
}

void ProbeSampler::sampling_locations(SampleLocType& locs) const
{
    locs.resize(m_npts);
    for (int i = 0; i < m_npts; ++i)
        for (int d = 0; d < AMREX_SPACEDIM; ++d) locs[i][d] = m_probes[i][d];
}

void ProbeSampler::sampling_orientations(SampleLocType& oris) const
{
    oris.resize(m_npts);
    for (int i = 0; i < m_npts; ++i)
        for (int d = 0; d < AMREX_SPACEDIM; ++d) oris[i][d] = m_oris[i][d];
}

#ifdef AMR_WIND_USE_NETCDF
void ProbeSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
{
    grp.put_attr("sampling_type", identifier());
}
#else
void ProbeSampler::define_netcdf_metadata(const ncutils::NCGroup&) const {}
#endif

} // namespace sampling
} // namespace amr_wind
