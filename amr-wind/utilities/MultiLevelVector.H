#ifndef MULTILEVELVECTOR_H
#define MULTILEVELVECTOR_H

#include "AMReX_Real.H"
#include "AMReX_Geometry.H"
#include "AMReX_Vector.H"
#include "AMReX_GpuContainers.H"

namespace amr_wind {

/** Multi-level vector
 *  \ingroup utilities
 *
 *  This class enables multi-level vectors. The primary vector is the
 *  size of the number of levels and it contains DeviceVectors size
 *  for the number of points in a given direction of the domain.
 */
class MultiLevelVector
{
public:
    MultiLevelVector() = default;

    void resize(const int axis, const amrex::Vector<amrex::Geometry>& geom);

    amrex::Gpu::DeviceVector<amrex::Real>& data(const int lev)
    {
        return m_data[lev];
    };

private:
    amrex::Vector<amrex::Real> m_dx;
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real>> m_data;
};
} // namespace amr_wind
#endif
