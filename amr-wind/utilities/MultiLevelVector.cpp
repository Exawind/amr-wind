#include "amr-wind/utilities/MultiLevelVector.H"

namespace amr_wind {

void MultiLevelVector::resize(
    const int axis, const amrex::Vector<amrex::Geometry>& geom)
{
    const auto nlevels = geom.size();
    m_data_h.resize(nlevels);
    m_data_d.resize(nlevels);
    for (int lev = 0; lev < nlevels; ++lev) {

        const int nc = geom[lev].Domain().length()[axis] +
                       ((m_floc == FieldLoc::CELL) ? 0 : 1);
        m_data_h[lev].resize(nc);
        m_data_h[lev].assign(nc, 0);
        m_data_d[lev].resize(nc);
    }
    copy_host_to_device();
};

void MultiLevelVector::copy_host_to_device()
{
    for (int lev = 0; lev < m_data_h.size(); ++lev) {
        amrex::Gpu::copyAsync(
            amrex::Gpu::hostToDevice, m_data_h[lev].begin(),
            m_data_h[lev].end(), m_data_d[lev].begin());
    }
}

} // namespace amr_wind
