#ifndef VOF_BCOP_H
#define VOF_BCOP_H

#include "amr-wind/equation_systems/vof/vof.H"
#include "amr-wind/equation_systems/BCOps.H"
#include "amr-wind/equation_systems/vof/vof_bciface.H"

namespace amr_wind::pde {

/** BC operations for VOF
 *  \ingroup vof
 */
template <>
struct BCOp<VOF>
{
    BCOp(PDEFields& fields, const SimTime& time)
        : m_fields(fields), m_time(time)
    {}

    void init_bcs()
    {
        BCVOF bc(m_fields.field);
        bc(VOF::default_bc_value);

        // No source terms currently exist for VOF, so this is inconsequential
        BCSrcTerm bc_src(m_fields.src_term);
        bc_src();
    }

    void apply_bcs(const FieldState rho_state)
    {
        amrex::IntVect ng_diff(1);
        auto& field = m_fields.field;
        const auto bc_time = (rho_state == FieldState::Old)
                                 ? m_time.current_time()
                                 : m_time.new_time();
        if ((rho_state != FieldState::Old && rho_state != FieldState::New)) {
            amrex::Abort(
                "vof_bcop.H apply_bcs(): a state other than New or Old was "
                "used. The method of setting bc_time must be evaluated before "
                "using a different state for this routine.");
        }
        field.fillphysbc(bc_time, ng_diff);
        field.apply_bc_funcs(rho_state);
    }

    PDEFields& m_fields;
    const SimTime& m_time;
};

} // namespace amr_wind::pde

#endif /* VOF_BCOP_H */
