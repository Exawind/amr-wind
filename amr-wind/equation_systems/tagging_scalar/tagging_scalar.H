#ifndef TAGGINGSCALAR_H
#define TAGGINGSCALAR_H

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/PDE.H"
#include "amr-wind/equation_systems/tagging_scalar/TaggingScalarSource.H"
#include "amr-wind/equation_systems/SourceTerm.H"

namespace amr_wind {
namespace pde {

/**
 *  \defgroup tagging_eqn Tagging scalar transport equation
 *  %Tagging scalar transport equation
 *
 *  \ingroup eqsys
 */

/**
 *  \defgroup temp_src Tagging Scalar Source terms
 *
 *  Source terms used with tagging scalar transport equation
 *
 *  \ingroup tagging_eqn
 */

/** Characteristics of tagging scalar transport equation
 *  \ingroup tagging_eqn
 */
struct TaggingScalar : ScalarTransport
{
    using MLDiffOp = amrex::MLABecLaplacian;
    using SrcTerm = SourceTerm; // TaggingScalarSource;

    static std::string pde_name() { return "TaggingScalar"; }
    static std::string var_name() { return "tagging_scalar"; }

    static constexpr amrex::Real default_bc_value = 0.0;

    static constexpr int ndim = 1;
    static constexpr bool multiply_rho = true;
    static constexpr bool has_diffusion = true;
    static constexpr bool need_nph_state = true;
};

/** Effective scalar diffusivity update operator
 *  \ingroup temp_eqn
 */
template <>
struct TurbulenceOp<TaggingScalar>
{
    TurbulenceOp(turbulence::TurbulenceModel& tmodel, PDEFields& fields)
        : m_tmodel(tmodel), m_fields(fields)
    {}

    void operator()()
    {
        auto& mueff = m_fields.mueff;
        m_tmodel.update_scalar_diff(mueff, TaggingScalar::var_name());
    }

    turbulence::TurbulenceModel& m_tmodel;
    PDEFields& m_fields;
};

} // namespace pde
} // namespace amr_wind

#endif /* TAGGINGSCALAR_H */
