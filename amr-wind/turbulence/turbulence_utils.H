#ifndef TURBULENCE_UTILS_H
#define TURBULENCE_UTILS_H

#include <memory>

#include "amr-wind/CFDSim.H"
#include "amr-wind/transport_models/ConstTransport.H"
#include "amr-wind/turbulence/LaminarModel.H"
#include "amr-wind/turbulence/LES/Smagorinsky.H"

namespace amr_wind {
namespace turbulence {

std::unique_ptr<TurbulenceModel> create_turbulence_model(
    const std::string& name, CFDSim& sim)
{
    std::unique_ptr<TurbulenceModel> turb;

    if (name == "Laminar") {
        turb.reset(new Laminar<transport::ConstTransport>(sim));
    }
    else if (name == "Smagorinsky") {
        turb.reset(new Smagorinsky<transport::ConstTransport>(sim));
    }
    else {
        amrex::Abort("Invalid turbulence model requested: " + name);
    }

    AMREX_ALWAYS_ASSERT(static_cast<bool>(turb));
    return turb;
}

}
} // namespace amr_wind

#endif /* TURBULENCE_UTILS_H */
