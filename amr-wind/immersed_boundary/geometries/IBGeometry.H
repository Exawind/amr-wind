#ifndef IBGEOMETRY_H
#define IBGEOMETRY_H

#include "amr-wind/immersed_boundary/IBtypes.H"
#include "amr-wind/core/vs/vector_space.H"

#include <string>

namespace amr_wind {
namespace ib {

/** Base data representation of an IB geometry.
 *
 *  \ingroup ib
 */
struct IBGeometryBaseData
{
    //! Number of points (markers) along the IB surface
    int num_pts;

    //! Body velocity at the boundary node
    VecList vel_body;

    //! Length of the IB boundary "panel"
    RealList dl;

    //! Total integrated lift over the wing
    amrex::Real lift{0.0};

    //! Total integrated drag over the wing
    amrex::Real drag{0.0};
};

struct IBGeometryType : public IBType
{};

} // namespace ib
} // namespace amr_wind

#endif /* IBGEOMETRY_H */
