#ifndef REAVERAGING_H
#define REAVERAGING_H

#include "amr-wind/utilities/averaging/TimeAveraging.H"

namespace amr_wind {
namespace averaging {

/** Reynolds averaging to compute mean field from a given CFD field
 */
class ReAveraging : public FieldTimeAverage::Register<ReAveraging>
{
public:
    static std::string identifier() { return "ReAveraging"; }

    ReAveraging(CFDSim&, const std::string& fname);

    /** Update field averaging at a given timestep
     *
     *  \param time SimTime instance
     *  \param filter_width Time-averaging window specified by user
     *
     *  \param elapsed_time Time elapsed since averaging was initiated
     */
    virtual void
    operator()(const SimTime&, const amrex::Real, const amrex::Real) override;

    virtual const std::string& average_field_name() override;

private:
    //! Generate the averaged field name based on the field name
    static std::string avg_name(const std::string& fname)
    {
        return fname + "_mean";
    }

    //! Fluctuating field
    const Field& m_field;

    //! Reynolds averaged field
    Field& m_average;
};

} // namespace averaging
} // namespace amr_wind

#endif /* REAVERAGING_H */
