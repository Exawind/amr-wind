#include "PDETraits.H"
#include "SchemeTraits.H"
#include "PDEHelpers.H"
#include "PDE.H"

namespace amr_wind {
namespace pde {

template class PDESystem<Temperature, fvm::Godunov>;
template class PDESystem<Temperature, fvm::MOL>;

} // namespace pde
} // namespace amr_wind
