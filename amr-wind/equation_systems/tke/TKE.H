#ifndef TKE_H
#define TKE_H

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/PDE.H"
#include "amr-wind/equation_systems/tke/TKESource.H"

namespace amr_wind {
namespace pde {

struct TKE : ScalarTransport
{
    using SrcTerm = TKESource;

    static std::string pde_name() { return "TKE"; }
    static std::string var_name() { return "tke"; }
};

} // namespace pde
} // namespace amr_wind

#endif /* TKE_H */
