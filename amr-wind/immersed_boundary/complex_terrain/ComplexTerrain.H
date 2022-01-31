#ifndef COMPLEX_TERRAIN_H
#define COMPLEX_TERRAIN_H

#include "amr-wind/immersed_boundary/IBTypes.H"
#include "amr-wind/core/vs/vector_space.H"

#include <string>

namespace amr_wind {
namespace ib {

/** Base data representation of an IB Complex Terrain geometry.
 *
 *  \ingroup ib
 */
struct ComplexTerrainBaseData
{

    bool is_mms{false};

    bool has_wall_model{false};

    amrex::Vector<amrex::Real> vel_bc{{0.0, 0.0, 0.0}};

    //! Total integrated drag forces on the immersed body
    amrex::Vector<amrex::Real> frc{{0.0, 0.0, 0.0}};
};

struct ComplexTerrainType : public IBType
{};

} // namespace ib
} // namespace amr_wind

#endif /* COMPLEX_TERRAIN_H */
