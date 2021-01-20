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
    using MLDiffOp = amrex::MLABecLaplacian;
    using SrcTerm = TKESource;

    static std::string pde_name() { return "TKE"; }
    static std::string var_name() { return "tke"; }

    static constexpr amrex::Real default_bc_value = 0.0;

    static constexpr int ndim = 1;
    static constexpr bool multiply_rho = true;
    static constexpr bool has_diffusion = true;
    static constexpr bool need_nph_state = true;
};

/** Implementation of Turbulence operator for TKE equation
 *  \ingroup tke_eqn
 */
template <>
struct TurbulenceOp<TKE>
{
    TurbulenceOp(turbulence::TurbulenceModel& tmodel, PDEFields& fields)
        : m_tmodel(tmodel), m_fields(fields)
    {}

    void operator()()
    {
        auto& mueff = m_fields.mueff;
        m_tmodel.update_scalar_diff(mueff, "tke");
    }

    turbulence::TurbulenceModel& m_tmodel;
    PDEFields& m_fields;
};

/** Implementation of post-processing actions for TKE equation
 *  \ingroup tke_eqn
 */
template <>
struct PostSolveOp<TKE>
{
    PostSolveOp(CFDSim& sim, PDEFields& fields) : m_sim(sim), m_fields(fields)
    {}

    void operator()(const amrex::Real time)
    {
        amrex::Print() << "Min/Max before lower bound op on TKE is "
                  << m_fields.field(0).min(0) << " - "
                  << m_fields.field(0).max(0) << std::endl;

        auto& field = m_fields.field;
        const auto& repo = field.repo();
        const int nlevels = repo.num_active_levels();
        for (int lev = 0; lev < nlevels; ++lev) {

            for (amrex::MFIter mfi(field(lev)); mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                const auto& field_arr = field(lev).array(mfi);

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        field_arr(i, j, k) = (field_arr(i, j, k) < 0.0)
                                                 ? clip_value
                                                 : field_arr(i, j, k);
                    });
            }
        }
        m_fields.field.fillpatch(time);
    }

    CFDSim& m_sim;
    PDEFields& m_fields;
    amrex::Real clip_value{1.0e-8};
};

} // namespace pde
} // namespace amr_wind

#endif /* TKE_H */
