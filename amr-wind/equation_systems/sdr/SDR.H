#ifndef SDR_H
#define SDR_H

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/PDE.H"
#include "amr-wind/equation_systems/sdr/SDRSource.H"

namespace amr_wind {
namespace pde {

/**
 *  \defgroup sdr_eqn Specific Dissipation Rate (SDR) transport equation
 *  %SDR transport equation
 *
 *  \ingroup eqsys
 */

/**
 *  \defgroup sdr_src SDR source terms
 *
 *  Source terms used with SDR transport equation
 *
 *  \ingroup sdr_eqn
 */

/** Characteristics of SDR transport equation
 *  \ingroup sdr_eqn
 */
struct SDR : ScalarTransport
{
    using SrcTerm = SDRSource;

    static std::string pde_name() { return "SDR"; }
    static std::string var_name() { return "sdr"; }
};

} // namespace pde
} // namespace amr_wind

#endif /* SDR_H */
