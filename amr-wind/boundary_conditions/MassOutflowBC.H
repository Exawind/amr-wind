#ifndef MASSINFLOWBC_H
#define MASSINFLOWBC_H

#include "amr-wind/core/FieldBCOps.H"
#include "amr-wind/core/FieldRepo.H"

#include "AMReX_Orientation.H"

namespace amr_wind {

/** Applies a dirichlet BC after having calculated du/dt+U*du/dx=0
 *  \ingroup field_bc
 *
 *  Used with the Dirichlet linear solver BC to apply a
 *  on-the-fly calculated dirichlet at the domain boundary.
 */
class MassOutflowBC : public FieldBCIface
{
public:
    MassOutflowBC(Field& field, amrex::Orientation ori);

    void operator()(Field& field, const FieldState /*rho_state*/) override;

private:
    Field& m_field;

    amrex::Orientation m_ori;
};

void MassOutflowBC::operator()(Field& /*field*/, const FieldState /*rho_state*/)
{}

} // namespace amr_wind

#endif /* MASSINFLOWBC_H */
