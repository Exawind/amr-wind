#ifndef LEVELSET_H
#define LEVELSET_H

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/PDE.H"

namespace amr_wind {
namespace pde {

struct Levelset : ScalarTransport
{

    static std::string pde_name() { return "Levelset"; }
    static std::string var_name() { return "levelset"; }

    static constexpr bool multiply_rho = false;
    static constexpr bool has_diffusion = false;

    static constexpr amrex::Real default_bc_value = 0.0;
};

} // namespace pde
} // namespace amr_wind

#endif /* LEVELSET_H */
