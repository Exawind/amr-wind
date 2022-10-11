#ifndef VOF_H
#define VOF_H

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/PDE.H"
#include "amr-wind/equation_systems/SourceTerm.H"

namespace amr_wind::pde {

struct VOF : ScalarTransport
{
    using MLDiffOp = amrex::MLABecLaplacian;
    using SrcTerm = SourceTerm;

    static std::string pde_name() { return "VOF"; }
    static std::string var_name() { return "vof"; }

    static constexpr int ndim = 1;
    static constexpr bool multiply_rho = false;
    static constexpr bool has_diffusion = false;
    static constexpr bool need_nph_state = true;

    static constexpr amrex::Real default_bc_value = 0.0;
};

} // namespace amr_wind::pde

#endif /* VOF_H */
