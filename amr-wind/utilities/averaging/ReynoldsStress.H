#ifndef REYNOLDSSTRESS_H
#define REYNOLDSSTRESS_H

#include "amr-wind/utilities/averaging/TimeAveraging.H"

namespace amr_wind {
namespace averaging {

/** Compute correlation from 2 CFD fields
 *  The formula used to compute the correlations between two fields:
 *  A = <A> + a and B = <B> + b is <ab> = <AB>-<A><B>
 *  where A and B are the mean values and a and b are the fluctuations
 *
 */
class ReynoldsStress : public FieldTimeAverage::Register<ReynoldsStress>
{
public:
    static std::string identifier() { return "ReynoldsStress"; }

    ReynoldsStress(CFDSim&, const std::string& fname);

    virtual void operator()(const SimTime&, const amrex::Real) override;

    virtual const std::string& average_field_name() override;

private:
    //! Generate the averaged field name based on the field name
    static std::string avg_name(const std::string& fname)
    {
        return "ReynoldsStress";
    }

    //! Fluctuating field
    const Field& m_field;

    //! Reynolds averaged field
    const Field& m_average;

    //! The stresses <AB>
    Field& m_stress;

    //! The reynolds stresses <ab>=<AB> - <A><B>
    Field& m_re_stress;
};

} // namespace averaging
} // namespace amr_wind

#endif /* REYNOLDSSTRESS_H */
