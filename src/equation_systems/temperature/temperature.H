#ifndef TEMPERATURE_H
#define TEMPERATURE_H

#include "PDETraits.H"
#include "SchemeTraits.H"
#include "PDEHelpers.H"
#include "PDE.H"

namespace amr_wind {
namespace pde {

struct Temperature : ScalarTransport
{
    static std::string pde_name() { return "Temperature"; }
    static std::string var_name() { return "temperature"; }
};

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
