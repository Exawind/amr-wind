#ifndef TIMEAVERAGING_H
#define TIMEAVERAGING_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/utilities/PostProcessing.H"

#include "AMReX_Vector.H"

#include <limits>
#include <memory>
#include <string>
#include <map>

namespace amr_wind {

class CFDSim;
class SimTime;
class Field;

namespace averaging {

/** Abstract class for time-averaging of CFD fields.
 *
 *  \ingroup utilities
 */
class FieldTimeAverage
    : public Factory<FieldTimeAverage, CFDSim&, const std::string&>
{
public:
    static std::string base_identifier() { return "FieldTimeAverage"; }

    virtual ~FieldTimeAverage() = default;

    /** Update field averaging at a given timestep
     *
     *  \param time SimTime instance
     *  \param filter_width Time-averaging window specified by user
     *
     *  \param elapsed_time Time elapsed since averaging was initiated
     */
    virtual void
    operator()(const SimTime&, const amrex::Real, const amrex::Real) = 0;

    virtual const std::string& average_field_name() = 0;
};

/** A collection of time-averaged quantities
 */
class TimeAveraging : public PostProcessBase::Register<TimeAveraging>
{
public:
    static std::string identifier() { return "TimeAveraging"; }

    TimeAveraging(CFDSim&, const std::string&);

    virtual ~TimeAveraging();

    //! Perform actions before mesh is created
    void pre_init_actions() override;

    void initialize() override;

    void post_advance_work() override;

    void post_regrid_actions() override {}

    const std::string& add_averaging(
        const std::string& field_name,
        const std::string& avg_type = "ReAveraging");

private:
    CFDSim& m_sim;

    const std::string m_label;

    //! Vector containing all registered field averagers
    amrex::Vector<std::unique_ptr<FieldTimeAverage>> m_averages;

    //! Fields registered so far to avoid duplication
    std::map<std::string, FieldTimeAverage*> m_registered;

    //! Time to start averaging the fields
    amrex::Real m_start_time{0.0};

    //! Time to stop averaging the fields
    amrex::Real m_stop_time{std::numeric_limits<amrex::Real>::max()};

    //! Time averaging window (in seconds)
    amrex::Real m_filter{0.0};
};

} // namespace averaging
} // namespace amr_wind

#endif /* TIMEAVERAGING_H */
