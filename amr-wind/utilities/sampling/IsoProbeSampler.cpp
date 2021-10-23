#include "amr-wind/utilities/sampling/IsoProbeSampler.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace sampling {

IsoProbeSampler::IsoProbeSampler(const CFDSim&) {}

IsoProbeSampler::~IsoProbeSampler() = default;

void IsoProbeSampler::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);
    std::string pfile("isoprobe_locations.txt");
    pp.query("isoprobe_location_file", pfile);
    std::string ofile("isoprobe_orientations.txt");
    pp.query("isoprobe_orientation_file", ofile);

    // Initial positions
    std::ifstream ifh(pfile, std::ios::in);
    if (!ifh.good())
        amrex::Abort("Cannot find isoprobe location file: " + pfile);

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
        amrex::Abort("Cannot find isoprobe orentation file: " + ofile);

    int m_opts;
    ifo >> m_opts;
    if (m_opts != m_npts)
        amrex::Abort("Isoprobe location number conflicts with orientations");

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

void IsoProbeSampler::sampling_locations(SampleLocType& locs) const
{
    locs.resize(m_npts);
    for (int i = 0; i < m_npts; ++i)
        for (int d = 0; d < AMREX_SPACEDIM; ++d) locs[i][d] = m_probes[i][d];
}

void IsoProbeSampler::sampling_orientations(SampleLocType& oris) const
{
    oris.resize(m_npts);
    for (int i = 0; i < m_npts; ++i)
        for (int d = 0; d < AMREX_SPACEDIM; ++d) oris[i][d] = m_oris[i][d];
}

#ifdef AMR_WIND_USE_NETCDF
void IsoProbeSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
{
    grp.put_attr("sampling_type", identifier());
}
#else
void IsoProbeSampler::define_netcdf_metadata(const ncutils::NCGroup&) const {}
#endif

} // namespace sampling
} // namespace amr_wind
