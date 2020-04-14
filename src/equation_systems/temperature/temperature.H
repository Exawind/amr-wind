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


} // namespace pde
} // namespace amr_wind


#endif /* TEMPERATURE_H */
