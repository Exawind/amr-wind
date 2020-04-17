#include "temperature/temperature.H"
#include "AdvOp_Godunov.H"
#include "AdvOp_MOL.H"

namespace amr_wind {
namespace pde {

template class PDESystem<Temperature, fvm::Godunov>;
template class PDESystem<Temperature, fvm::MOL>;

} // namespace pde
} // namespace amr_wind
