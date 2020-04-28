#ifndef SDR_H
#define SDR_H

#include "PDETraits.H"
#include "SchemeTraits.H"
#include "PDEHelpers.H"
#include "PDE.H"
#include "sdr/SDRSource.H"

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
