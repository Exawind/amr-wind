#ifndef MULTILEVELVECTOR_H
#define MULTILEVELVECTOR_H

#include "AMReX_Geometry.H"
#include "AMReX_Vector.H"
#include "AMReX_GpuContainers.H"

namespace amr_wind {

/** Multi-level vector
 *  \ingroup utilities
 *
 *  This class enables multi-level vectors. The primary vector is the
 *  size of the number of levels and it contains host and device
 *  vectors sized for the number of points in a given direction of the
 *  domain.
 */
class MultiLevelVector
{
public:
    MultiLevelVector() = default;

    void resize(const int axis, const amrex::Vector<amrex::Geometry>& geom);

    amrex::Vector<amrex::Real>& host_data(const int lev)
    {
        return m_data_h[lev];
    };

    amrex::Gpu::DeviceVector<amrex::Real>& device_data(const int lev)
    {
        return m_data_d[lev];
    };

    void copy_host_to_device();

    int size() const { return static_cast<int>(m_data_h.size()); }

    int ncells(const int lev) const
    {
        return static_cast<int>(m_data_h[lev].size());
    }

private:
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real>> m_data_d;
    amrex::Vector<amrex::Vector<amrex::Real>> m_data_h;
};
} // namespace amr_wind
#endif
