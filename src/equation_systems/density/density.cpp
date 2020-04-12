#include "PDETraits.H"
#include "SchemeTraits.H"
#include "PDEHelpers.H"
#include "PDE.H"

namespace amr_wind {
namespace pde {

template class PDESystem<Density, fvm::Godunov>;
template class PDESystem<Density, fvm::MOL>;

} // namespace pde
} // namespace amr_wind
