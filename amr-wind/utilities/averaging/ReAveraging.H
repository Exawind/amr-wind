#ifndef REAVERAGING_H
#define REAVERAGING_H

#include "amr-wind/utilities/averaging/TimeAveraging.H"

namespace amr_wind::averaging {

/** Reynolds averaging to compute mean field from a given CFD field
 */
class ReAveraging : public FieldTimeAverage::Register<ReAveraging>
{
public:
    static std::string identifier() { return "ReAveraging"; }

    ReAveraging(CFDSim& /*sim*/, const std::string& avgname, const std::string& fname);

    /** Update field averaging at a given timestep
     *
     *  \param time SimTime instance
     *  \param filter_width Time-averaging window specified by user
     *
     *  \param elapsed_time Time elapsed since averaging was initiated
     */
    void operator()(
        const SimTime& /*time*/,
        const amrex::Real /*filter_width*/,
        const int /*avg_frequency*/,
        const amrex::Real /*elapsed_time*/) override;

    const std::string& average_field_name() override;

private:
    //! Generate the averaged field name based on the field name
    static std::string avg_name(const std::string& fname, 
                                const std::string& avgname)
    {
        return fname + "_" + "mean" + "_" + avgname;
      //return fname + "_" + "mean";
    }

    //! Fluctuating field
    const Field& m_field;

    //! Reynolds averaged field
    Field& m_average;

    //! Averaging instance name
    const std::string& m_avgname;
};

} // namespace amr_wind::averaging

#endif /* REAVERAGING_H */
