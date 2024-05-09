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

private:
    // Functions called within public functions
    void update_gradp();

    CFDSim* m_sim_ptr;
};

} // namespace amr_wind

#endif /* OVERSETOPS_H */