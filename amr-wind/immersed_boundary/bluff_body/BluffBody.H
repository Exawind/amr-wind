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
