#ifndef DENSITY_H
#define DENSITY_H

#include "PDETraits.H"
#include "SchemeTraits.H"
#include "PDEHelpers.H"
#include "PDE.H"

namespace amr_wind {
namespace pde {

struct Density : ScalarTransport
{
    static std::string pde_name() { return "Density"; }
    static std::string var_name() { return "density"; }

    static constexpr bool multiply_rho = false;
};

} // namespace pde
} // namespace amr_wind

#endif /* DENSITY_H */
