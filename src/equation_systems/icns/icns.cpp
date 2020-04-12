#include "PDETraits.H"
#include "SchemeTraits.H"
#include "PDEHelpers.H"
#include "PDE.H"
#include "icns/icns_ops.H"

namespace amr_wind {
namespace pde {

template class PDESystem<ICNS, ::amr_wind::fvm::Godunov>;
template class PDESystem<ICNS, ::amr_wind::fvm::MOL>;

} // namespace pde
} // namespace amr_wind
