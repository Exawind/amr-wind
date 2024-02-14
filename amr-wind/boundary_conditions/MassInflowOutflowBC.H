#ifndef MASSINFLOWOUTFLOWBC_H
#define MASSINFLOWOUTFLOWBC_H

#include "amr-wind/core/FieldBCOps.H"
#include "amr-wind/core/FieldRepo.H"

#include "AMReX_Orientation.H"

namespace amr_wind {

/** Can handle both inflow and outflow cells in the same boundary
 *  \ingroup field_bc
 *
 *
 */
class MassInflowOutflowBC : public FieldBCIface
{
public:
    MassInflowOutflowBC(Field& field, amrex::Orientation ori);

    void operator()(Field& field, const FieldState /*rho_state*/) override;

private:
    Field& m_field;

    amrex::Orientation m_ori;
};

} // namespace amr_wind


#endif /* MASSINFLOWOUTFLOWBC_H */

