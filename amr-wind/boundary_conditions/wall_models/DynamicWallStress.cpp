#include "amr-wind/boundary_conditions/wall_models/DynamicWallStress.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/diffusion/diffusion.H"

#include <cmath>

#include "AMReX_ParmParse.H"
#include "AMReX_Print.H"
#include "AMReX_ParallelDescriptor.H"

namespace amr_wind {

DynamicWallStress::DynamicWallStress(CFDSim& sim)
    : m_sim(sim), m_mesh(m_sim.mesh())
{
    // Here shouldbe all the wave-related stuff
}
} // namespace amr_wind
