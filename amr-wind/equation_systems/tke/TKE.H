#ifndef TKE_H
#define TKE_H

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/PDE.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/equation_systems/tke/TKESource.H"

namespace amr_wind {
namespace pde {

/**
 *  \defgroup tke_eqn Turbulent Kinetic Energy (TKE) transport equation
 *  %TKE transport equation
 *
 *  \ingroup eqsys
 */

/**
 *  \defgroup tke_src TKE Source terms
 *
 *  Source terms used with TKE transport equation
 *
 *  \ingroup tke_eqn
 */

/** Characteristics of TKE transport equation
 *
 *  \ingroup tke_eqn
 */
struct TKE : ScalarTransport
{
    using SrcTerm = TKESource;

    static std::string pde_name() { return "TKE"; }
    static std::string var_name() { return "tke"; }
};

/** Implementation of Turbulence operator for TKE equation
 *  \ingroup tke_eqn
 */
template<>
struct TurbulenceOp<TKE>
{
    TurbulenceOp(turbulence::TurbulenceModel& tmodel,
                 PDEFields& fields)
        : m_tmodel(tmodel), m_fields(fields)
    {}

    void operator()()
    {
        auto& tvisc = m_tmodel.mu_turb();
        auto& mueff = m_fields.mueff;
        mueff.setVal(0.0);
        field_ops::saxpy(
            mueff, 2.0, tvisc, 0, 0, mueff.num_comp(), mueff.num_grow());
    }

    turbulence::TurbulenceModel& m_tmodel;
    PDEFields& m_fields;
};

/** Implementation of post-processing actions for TKE equation
 *  \ingroup tke_eqn
 */
template<>
struct PostSolveOp<TKE>
{
    PostSolveOp(PDEFields& fields) : m_fields(fields) {}

    void operator()(const amrex::Real time)
    {
        field_ops::lower_bound(m_fields.field, clip_value);
        m_fields.field.fillpatch(time);
    }

    PDEFields& m_fields;
    amrex::Real clip_value{1.0e-15};
};

} // namespace pde
} // namespace amr_wind

#endif /* TKE_H */
