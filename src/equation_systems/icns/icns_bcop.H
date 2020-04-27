#ifndef ICNS_BCOP_H
#define ICNS_BCOP_H

#include "icns/icns.H"
#include "BCOps.H"

namespace amr_wind {
namespace pde {

/** BC operations for incompressible Navier-Stokes
 */
template<>
struct BCOp<ICNS>
{
    BCOp(PDEFields& fields, const SimTime& time)
        : m_fields(fields), m_time(time)
    {}

    /** Initialize the BCs for the PDE variable as well as its source term
     */
    void init_bcs(const BCIface::IncfloBC& ibc_type)
    {
        BCVelocity bc(m_fields.field, ibc_type);
        bc();

        auto& density = m_fields.repo.get_field("density");
        const amrex::Real density_default = 1.0;
        BCScalar bc_den(density, ibc_type);
        bc_den(density_default);

        auto& pressure = m_fields.repo.get_field("p");
        BCPressure bc_pres(pressure, ibc_type);
        bc_pres();

        // Used for fillpatch operation on the source term
        BCSrcTerm bc_src(m_fields.src_term, ibc_type);
        bc_src();
    }

    /** Apply boundary conditions before a linear solve
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

}
}

#endif /* ICNS_BCOP_H */
