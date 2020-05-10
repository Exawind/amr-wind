#ifndef FIXEDGRADIENTBC_H
#define FIXEDGRADIENTBC_H

#include "FieldBCOps.H"
#include "FieldRepo.H"

#include "AMReX_Orientation.H"

namespace amr_wind {

class FixedGradientBC : public FieldBCIface
{
public:
    FixedGradientBC(Field& field, amrex::Orientation ori);

    void operator()(Field& field, const FieldState) override;

private:
    Field& m_field;

    amrex::Orientation m_ori;
};

}

#endif /* FIXEDGRADIENTBC_H */
