#include "sdr/SDR.H"
#include "AdvOp_Godunov.H"
#include "AdvOp_MOL.H"
#include "BCOps.H"

namespace amr_wind {
namespace pde {

template class PDESystem<SDR, fvm::Godunov>;
template class PDESystem<SDR, fvm::MOL>;

} // namespace pde
} // namespace amr_wind
