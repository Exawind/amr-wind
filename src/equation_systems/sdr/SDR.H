#ifndef SDR_H
#define SDR_H

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/PDE.H"
#include "amr-wind/equation_systems/sdr/SDRSource.H"

namespace amr_wind {
namespace pde {

struct SDR : ScalarTransport
{
    using SrcTerm = SDRSource;

    static std::string pde_name() { return "SDR"; }
    static std::string var_name() { return "sdr"; }
};

} // namespace pde
} // namespace amr_wind

#endif /* SDR_H */
