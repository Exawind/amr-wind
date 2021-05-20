#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "amr-wind/immersed_boundary/IBTypes.H"
#include "amr-wind/core/vs/vector_space.H"

#include <string>

namespace amr_wind {
namespace ib {

/** Base data representation of an IB geometry.
 *
 *  \ingroup ib
 */
struct BluffBodyBaseData
{
    //! Number of points (markers) along the IB surface
    int num_pts;

    //! Body velocity at the boundary node
    VecList body_vel;

    //! Forcing volume in theory should be equal to dx[0]*dx[1]*dx[2]
    RealList dVl;

    //! Total integrated lift over the wing
    amrex::Real lift{0.0};

    //! Total integrated drag over the wing
    amrex::Real drag{0.0};
};

struct BluffBodyType : public IBType
{};

} // namespace ib
} // namespace amr_wind

#endif /* GEOMETRY_H */
