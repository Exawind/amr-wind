#include "amr-wind/utilities/MultiLevelVector.H"

namespace amr_wind {
void MultiLevelVector::resize(
    const int axis, const amrex::Vector<amrex::Geometry>& geom)
{
    const auto nlevels = geom.size();
    m_data.resize(nlevels);
    m_dx.resize(nlevels);
    for (int lev = 0; lev < nlevels; ++lev) {
        m_dx[lev] = geom[lev].CellSize()[axis];

        const int ncells = geom[lev].Domain().length()[axis];
        m_data[lev].resize(ncells);
        m_data[lev].assign(ncells, 0);
    }
};

} // namespace amr_wind
