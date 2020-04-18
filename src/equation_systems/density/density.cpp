#include "density/density.H"
#include "AdvOp_Godunov.H"
#include "AdvOp_MOL.H"
#include "density/density_ops.H"

namespace amr_wind {
namespace pde {

template class PDESystem<Density, fvm::Godunov>;
template class PDESystem<Density, fvm::MOL>;

} // namespace pde
} // namespace amr_wind
