#include "amr-wind/equation_systems/passive_tracer/passive_tracer.H"
#include "amr-wind/equation_systems/AdvOp_Godunov.H"
#include "amr-wind/equation_systems/AdvOp_MOL.H"
#include "amr-wind/equation_systems/BCOps.H"

namespace amr_wind {
namespace pde {

template class PDESystem<PassiveTracer, fvm::Godunov>;
template class PDESystem<PassiveTracer, fvm::MOL>;

} // namespace pde
} // namespace amr_wind
