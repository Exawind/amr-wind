#ifndef TKE_H
#define TKE_H

#include "PDETraits.H"
#include "SchemeTraits.H"
#include "PDEHelpers.H"
#include "PDE.H"
#include "tke/TKESource.H"

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
