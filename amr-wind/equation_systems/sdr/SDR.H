#ifndef SDR_H
#define SDR_H

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/PDE.H"
#include "amr-wind/core/field_ops.H"
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
    using MLDiffOp = amrex::MLABecLaplacian;
    using SrcTerm = SDRSource;

    static std::string pde_name() { return "SDR"; }
    static std::string var_name() { return "sdr"; }

    static constexpr amrex::Real default_bc_value = 0.0;

    static constexpr int ndim = 1;
    static constexpr bool multiply_rho = true;
    static constexpr bool has_diffusion = true;
    static constexpr bool need_nph_state = true;
};

template<>
struct TurbulenceOp<SDR>
{
    TurbulenceOp(turbulence::TurbulenceModel& tmodel,
                 PDEFields& fields)
        : m_tmodel(tmodel), m_fields(fields)
    {}

    void operator()()
    {
        auto& mueff = m_fields.mueff;
        m_tmodel.update_scalar_diff(mueff, "sdr");
    }

    turbulence::TurbulenceModel& m_tmodel;
    PDEFields& m_fields;
};

/** Implementation of post-processing actions for SDR equation
 *  \ingroup sdr_eqn
 */
template <>
struct PostSolveOp<SDR>
{
    PostSolveOp(PDEFields& fields) : m_fields(fields) {}
    
    void operator()(const amrex::Real time)
    {
        field_ops::lower_bound(m_fields.field, clip_value);
        m_fields.field.fillpatch(time);
    }

    PDEFields& m_fields;
    amrex::Real clip_value{1.0e-8};
};
    

} // namespace pde
} // namespace amr_wind

#endif /* SDR_H */
