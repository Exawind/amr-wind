#ifndef PassiveTracer_H
#define PassiveTracer_H

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/PDE.H"
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
struct PassiveTracer : ScalarTransport
{
    using MLDiffOp = amrex::MLABecLaplacian;
    using SrcTerm = SourceTerm;

    static std::string pde_name() { return "PassiveTracer"; }
    static std::string var_name() { return "passive_tracer"; }

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
struct TurbulenceOp<PassiveTracer>
{
    TurbulenceOp(turbulence::TurbulenceModel& tmodel, PDEFields& fields)
        : m_tmodel(tmodel), m_fields(fields)
    {}

    void operator()()
    {
        auto& mueff = m_fields.mueff;
        m_tmodel.update_scalar_diff(mueff, PassiveTracer::var_name());
    }

    turbulence::TurbulenceModel& m_tmodel;
    PDEFields& m_fields;
};

} // namespace pde
} // namespace amr_wind

#endif /* PassiveTracer_H */
