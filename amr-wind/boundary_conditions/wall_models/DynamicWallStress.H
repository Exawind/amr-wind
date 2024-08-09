#ifndef DYNAMIC_WALL_STRESS_H
#define DYNAMIC_WALL_STRESS_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldBCOps.H"

namespace amr_wind {

/** Wall-function computations for LES simulations
 *  \ingroup utilities
 *
 *  This class performs the necessary computations at the beginning of
 *  predictor/corrector steps. The actual BC population in ghost cells is
 *  performed by VelWallFunc BC interface class.
 */
class DynamicWallStress
{
public:
    explicit DynamicWallStress(CFDSim& sim);

    amrex::Real dyn_tau() const { return 0.; }

    ~DynamicWallStress() = default;

private:
    const CFDSim& m_sim;

    const amrex::AmrCore& m_mesh;

    int m_direction{2}; ///< Direction normal to wall, hardcoded to z
};

} // namespace amr_wind
#endif