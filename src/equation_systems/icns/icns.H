#ifndef ICNS_H
#define ICNS_H

#include "PDETraits.H"
#include "SchemeTraits.H"
#include "PDEHelpers.H"
#include "PDE.H"
#include "icns/MomentumSource.H"

namespace amr_wind {
namespace pde {

struct ICNS : VectorTransport
{
    using MLDiffOp = amrex::MLTensorOp;
    using SrcTerm = MomentumSource;

    static std::string pde_name() { return "ICNS"; }
    static std::string var_name() { return "velocity"; }

    static constexpr bool multiply_rho = false;
    static constexpr bool has_diffusion = true;

    // No n+1/2 state for velocity for now
    static constexpr bool need_nph_state = false;
};

} // namespace pde
} // namespace amr_wind

#endif /* ICNS_H */
