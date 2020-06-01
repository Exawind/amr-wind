#ifndef BCINTERFACE_H
#define BCINTERFACE_H

#include "amr-wind/core/Field.H"
#include "amr-wind/incflo_enums.H"

#include "AMReX_BC_TYPES.H"

namespace amr_wind {

/** Abstract class that defines the interface for processing BCs
 */
class BCIface
{
public:
    using IncfloBC = amrex::GpuArray<BC, AMREX_SPACEDIM * 2>;

    BCIface(Field& field);

    virtual ~BCIface() = default;

    //! Operator that performs init actions and syncs the BC data to device
    virtual void operator()(const amrex::Real value = 0.0);

protected:
    //! Setup AMReX mathematical BC types
    virtual void set_bcrec() = 0;

    //! Read input data for specific fields based on BC types
    virtual void read_values() = 0;

    //! Parse custom BCs specified for this field
    virtual void read_bctype();

    virtual void set_bcfuncs();

    //! Set default BC values for the field
    inline void set_default_value(const amrex::Real value);

    //! Set AMReX mathematical boundary types for the lower boundaries
    inline void set_bcrec_lo(
        int dir, amrex::BCType::mathematicalBndryTypes bcrec);

    //! Set AMReX mathematical boundary types for the upper boundaries
    inline void set_bcrec_hi(
        int dir, amrex::BCType::mathematicalBndryTypes bcrec);

    //! Field instance where BC is being set
    Field& m_field;
};

class BCVelocity : public BCIface
{
public:
    BCVelocity(Field& field) : BCIface(field) {}

protected:
    virtual void set_bcrec() override;
    virtual void read_values() override;
};

class BCScalar : public BCIface
{
public:
    BCScalar(Field& field) : BCIface(field) {}

protected:
    virtual void set_bcrec() override;
    virtual void read_values() override;
};

class BCPressure : public BCScalar
{
public:
    BCPressure(Field& field) : BCScalar(field) {}

protected:
    virtual void read_values() override;
};

class BCSrcTerm : public BCIface
{
public:
    BCSrcTerm(Field& field) : BCIface(field) {}

protected:
    virtual void set_bcrec() override;
    virtual void read_values() override
    {
        // No values specified for source terms. empty method to satisfy inheritance
    }
};

class BCFillPatchExtrap : public BCIface
{
public:
    BCFillPatchExtrap(
        Field& field,
        amrex::BCType::mathematicalBndryTypes bctype = amrex::BCType::hoextrap)
        : BCIface(field)
        , m_extrap_type(bctype)
    {}

protected:
    virtual void set_bcrec() override;
    virtual void read_values() override
    {
        // No values specified for source terms. empty method to satisfy
        // inheritance
    }

    amrex::BCType::mathematicalBndryTypes m_extrap_type;
};

} // namespace amr_wind

#endif /* BCINTERFACE_H */
