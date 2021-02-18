#include "amr-wind/utilities/averaging/TimeAveraging.H"
#include "amr-wind/utilities/averaging/ReAveraging.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace averaging {

TimeAveraging::TimeAveraging(CFDSim& sim, const std::string& label)
    : m_sim(sim), m_label(label)
{}

TimeAveraging::~TimeAveraging() = default;

void TimeAveraging::pre_init_actions()
{
    //! Different averaging types
    amrex::Vector<std::string> labels;
    {
        amrex::ParmParse pp(m_label);
        pp.getarr("labels", labels);
    }

    for (const auto& lbl : labels) {
        //! Fields to be averaged
        amrex::Vector<std::string> fnames;
        std::string avg_type;
        const std::string pp_key = m_label + "/" + lbl;

        amrex::ParmParse pp1(pp_key);
        pp1.getarr("fields", fnames);
        pp1.get("averaging_type", avg_type);
        pp1.get("averaging_window", m_filter);

        for (const auto& fname : fnames) {
            const std::string key = fname + "_" + avg_type;

            // Guard against multiple registrations of the field
            const auto found = m_registered.find(key);
            if (found != m_registered.end()) continue;

            // Create the averaging entity
            m_averages.emplace_back(
                FieldTimeAverage::create(avg_type, m_sim, fname));

            // Track fields that have an average
            m_registered.emplace(key, m_averages.back().get());
        }
    }
}

void TimeAveraging::initialize() {}

const std::string& TimeAveraging::add_averaging(
    const std::string& field_name, const std::string& avg_type)
{
    const std::string key = field_name + "_" + avg_type;
    const auto found = m_registered.find(key);

    // If the field was already registered then just return the field
    if (found != m_registered.end())
        return m_registered[key]->average_field_name();

    // Create and register new average
    m_averages.emplace_back(
        FieldTimeAverage::create(avg_type, m_sim, field_name));
    return m_averages.back()->average_field_name();
}

void TimeAveraging::post_advance_work()
{
    for (auto& avg : m_averages) {
        (*avg)(m_sim.time(), m_filter);
    }
}

} // namespace averaging
} // namespace amr_wind
