#include "amr-wind/equation_systems/tke/TKE.H"
#include "amr-wind/equation_systems/AdvOp_Godunov.H"
#include "amr-wind/equation_systems/AdvOp_MOL.H"
#include "amr-wind/equation_systems/BCOps.H"

namespace amr_wind {
namespace pde {

template class PDESystem<TKE, fvm::Godunov>;
template class PDESystem<TKE, fvm::MOL>;

} // namespace pde
} // namespace amr_wind
