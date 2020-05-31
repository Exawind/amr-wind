#include "amr-wind/equation_systems/temperature/temperature.H"
#include "amr-wind/equation_systems/AdvOp_Godunov.H"
#include "amr-wind/equation_systems/AdvOp_MOL.H"
#include "amr-wind/equation_systems/BCOps.H"

namespace amr_wind {
namespace pde {

template class PDESystem<Temperature, fvm::Godunov>;
template class PDESystem<Temperature, fvm::MOL>;

} // namespace pde
} // namespace amr_wind
