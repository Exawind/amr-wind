#ifndef TEMPERATURE_H
#define TEMPERATURE_H

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/PDE.H"
#include "amr-wind/equation_systems/temperature/TemperatureSource.H"

namespace amr_wind {
namespace pde {

/**
 *  \defgroup temp_eqn Temperature transport equation
 *  %Temperature transport equation
 *
 *  \ingroup eqsys
 */

/**
 *  \defgroup temp_src Temperature Source terms
 *
 *  Source terms used with temperature transport equation
 *
 *  \ingroup temp_eqn
 */

/** Characteristics of temperature transport equation
 *  \ingroup temp_eqn
 */
struct Temperature : ScalarTransport
{
    using SrcTerm = TemperatureSource;

    static std::string pde_name() { return "Temperature"; }
    static std::string var_name() { return "temperature"; }
};

/** Effective thermal diffusivity update operator
 *  \ingroup temp_eqn
 */
template<>
struct TurbulenceOp<Temperature>
{
    TurbulenceOp(turbulence::TurbulenceModel& tmodel,
                 PDEFields& fields)
        : m_tmodel(tmodel), m_fields(fields)
    {}

    void operator()()
    {
        m_tmodel.update_alphaeff(m_fields.mueff);
    }

    turbulence::TurbulenceModel& m_tmodel;
    PDEFields& m_fields;
};

} // namespace pde
} // namespace amr_wind


#endif /* TEMPERATURE_H */
