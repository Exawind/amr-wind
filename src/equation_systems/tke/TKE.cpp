#include "tke/TKE.H"
#include "AdvOp_Godunov.H"
#include "AdvOp_MOL.H"
#include "BCOps.H"

namespace amr_wind {
namespace pde {

template class PDESystem<TKE, fvm::Godunov>;
template class PDESystem<TKE, fvm::MOL>;

} // namespace pde
} // namespace amr_wind
