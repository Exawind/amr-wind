#ifndef BCOPS_H
#define BCOPS_H

#include "PDETraits.H"
#include "PDEOps.H"
#include "BCInterface.H"

namespace amr_wind {
namespace pde {

/** BC operations for all scalar transport equations
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
    void init_bcs(const BCIface::IncfloBC ibc_type)
    {
        BCScalar bc(m_fields.field, ibc_type);
        bc(PDE::default_bc_value);

        // Used for fillpatch operation on the source term
        BCScalar bc_src(m_fields.src_term, ibc_type);
        bc_src();
    }

    /**
     */
    void apply_bcs()
    {
        amrex::IntVect ng_diff(1);
        auto& field = m_fields.field;
        field.fillphysbc(m_time.new_time(), ng_diff);
    }

    PDEFields& m_fields;
    const SimTime& m_time;
};

} // namespace pde
} // namespace amr_wind

#endif /* BCOPS_H */
