#ifndef BCOPS_H
#define BCOPS_H

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/PDEOps.H"
#include "amr-wind/boundary_conditions/BCInterface.H"
#include "amr-wind/boundary_conditions/scalar_bcs.H"

namespace amr_wind::pde {

/** BC operations for all scalar transport equations
 *  \ingroup pdeop
 */
template <typename PDE>
struct BCOp<
    PDE,
    typename std::enable_if<std::is_base_of<ScalarTransport, PDE>::value>::type>
{
    BCOp(PDEFields& fields, const SimTime& time)
        : m_fields(fields), m_time(time)
    {}

    /** Initialize the BCs for the PDE variable as well as its source term
     */
    void init_bcs()
    {
        BCScalar bc(m_fields.field);
        bc(PDE::default_bc_value);
        const auto udfs = bc.get_dirichlet_udfs();
        scalar_bc::register_scalar_dirichlet(
            m_fields.field, m_fields.repo.mesh(), m_time, udfs);

        // Used for fillpatch operation on the source term
        BCSrcTerm bc_src(m_fields.src_term);
        bc_src();
    }

    /**
     */
    void apply_bcs(const FieldState rho_state)
    {
        amrex::IntVect ng_diff(1);
        auto& field = m_fields.field;
        field.fillphysbc(m_time.new_time(), ng_diff);
        field.apply_bc_funcs(rho_state);
    }

    PDEFields& m_fields;
    const SimTime& m_time;
};

} // namespace amr_wind::pde

#endif /* BCOPS_H */
