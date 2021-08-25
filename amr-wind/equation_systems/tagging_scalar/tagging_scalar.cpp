#include "amr-wind/equation_systems/tagging_scalar/tagging_scalar.H"
#include "amr-wind/equation_systems/AdvOp_Godunov.H"
#include "amr-wind/equation_systems/AdvOp_MOL.H"
#include "amr-wind/equation_systems/BCOps.H"

namespace amr_wind {
namespace pde {

template class PDESystem<TaggingScalar, fvm::Godunov>;
template class PDESystem<TaggingScalar, fvm::MOL>;

} // namespace pde
} // namespace amr_wind
