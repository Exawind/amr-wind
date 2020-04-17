#include "icns/icns.H"
#include "icns/icns_ops.H"
#include "icns/icns_diffusion.H"

namespace amr_wind {
namespace pde {

template class PDESystem<ICNS, ::amr_wind::fvm::Godunov>;
template class PDESystem<ICNS, ::amr_wind::fvm::MOL>;

} // namespace pde
} // namespace amr_wind
