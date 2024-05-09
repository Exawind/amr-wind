#ifndef OVERSETOPS_H
#define OVERSETOPS_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"

namespace amr_wind {

class OversetOps
{
public:
    OversetOps() = default;

    void initialize(CFDSim& sim);
    void pre_advance_actions();

    // Needs to be public for cuda
    void update_gradp();

private:
    CFDSim* m_sim_ptr;
};

} // namespace amr_wind

#endif /* OVERSETOPS_H */