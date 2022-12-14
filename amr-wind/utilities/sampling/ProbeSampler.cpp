#include "amr-wind/utilities/sampling/ProbeSampler.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::sampling {

ProbeSampler::ProbeSampler(const CFDSim& /*unused*/) {}

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

    ifh >> m_npts;
    ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    m_probes.resize(m_npts);
    for (int i = 0; i < m_npts; ++i) {
        ifh >> m_probes[i][0] >> m_probes[i][1] >> m_probes[i][2];
        ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
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
}
#else
void ProbeSampler::define_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
#endif

} // namespace amr_wind::sampling
