#include "amr-wind/utilities/MultiLevelVector.H"

namespace amr_wind {
void MultiLevelVector::resize(
    const int axis, const amrex::Vector<amrex::Geometry>& geom)
{
    const auto nlevels = geom.size();
    m_data_h.resize(nlevels);
    m_data_d.resize(nlevels);
    for (int lev = 0; lev < nlevels; ++lev) {

        const int ncells = geom[lev].Domain().length()[axis];
        m_data_h[lev].resize(ncells);
        m_data_h[lev].assign(ncells, 0);
        m_data_d[lev].resize(ncells);
    }
    sync_host_to_device();
};

void MultiLevelVector::sync_host_to_device()
{
    for (int lev = 0; lev < m_data_h.size(); ++lev) {
#ifdef AMREX_USE_GPU
        amrex::Gpu::htod_memcpy_async(
            m_data_d[lev].data(), m_data_h[lev].data(),
            sizeof(amrex::Real) * m_data_h[lev].size());
#else
        std::memcpy(
            m_data_d[lev].data(), m_data_h[lev].data(),
            sizeof(amrex::Real) * m_data_h[lev].size());
#endif
    }
}

} // namespace amr_wind
