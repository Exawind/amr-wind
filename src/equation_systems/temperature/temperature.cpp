#include "temperature/temperature.H"

namespace amr_wind {
namespace pde {

template class PDESystem<Temperature, fvm::Godunov>;
template class PDESystem<Temperature, fvm::MOL>;

} // namespace pde
} // namespace amr_wind
