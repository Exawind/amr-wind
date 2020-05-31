#include "amr-wind/equation_systems/icns/icns.H"
#include "amr-wind/equation_systems/icns/icns_ops.H"
#include "amr-wind/equation_systems/icns/icns_diffusion.H"
#include "amr-wind/equation_systems/icns/icns_bcop.H"

namespace amr_wind {
namespace pde {

template class PDESystem<ICNS, ::amr_wind::fvm::Godunov>;
template class PDESystem<ICNS, ::amr_wind::fvm::MOL>;

} // namespace pde
} // namespace amr_wind
