#ifndef ACTUATORWING_H
#define ACTUATORWING_H

#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/core/vs/vector_space.H"

#include <string>

namespace amr_wind {
namespace actuator {

/** Base data representation of a fixed-wing.
 *
 *  \ingroup actuator
 */
struct WingBaseData
{
    //! Number of points along the wing
    int num_pts;

    //! Starting coordinate of the wing
    vs::Vector start;

    //! Ending coordinate of the wing
    vs::Vector end;

    //! Normal to the wing
    vs::Vector normal{0.0, 0.0, 1.0};

    //! Gaussian smearing factor input by user
    vs::Vector eps_inp;

    //! Pitch angle for the wing
    amrex::Real pitch{0.0};
};

} // namespace actuator
} // namespace amr_wind

#endif /* ACTUATORWING_H */
