#include "amr-wind/equation_systems/sdr/SDR.H"
#include "amr-wind/equation_systems/AdvOp_Godunov.H"
#include "amr-wind/equation_systems/AdvOp_MOL.H"
#include "amr-wind/equation_systems/BCOps.H"

namespace amr_wind {
namespace pde {

template class PDESystem<SDR, fvm::Godunov>;
template class PDESystem<SDR, fvm::MOL>;

} // namespace pde
} // namespace amr_wind
