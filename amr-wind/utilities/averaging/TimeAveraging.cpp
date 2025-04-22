#include <utility>

#include "amr-wind/utilities/averaging/TimeAveraging.H"
#include "amr-wind/utilities/averaging/ReAveraging.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::averaging {

TimeAveraging::TimeAveraging(CFDSim& sim, std::string label)
    : m_sim(sim), m_label(std::move(label))
{}

TimeAveraging::~TimeAveraging() = default;

void TimeAveraging::pre_init_actions()
{
    //! Different averaging types
    amrex::Vector<std::string> labels;
    {
        amrex::ParmParse pp(m_label);
        pp.getarr("labels", labels);
        ioutils::assert_with_message(
            ioutils::all_distinct(labels),
            "Duplicates in " + m_label + ".labels");
        pp.query("averaging_start_time", m_start_time);
        pp.query("averaging_stop_time", m_stop_time);
        pp.query("averaging_frequency", m_frequency);
        pp.get("averaging_window", m_filter);
    }

    std::cout << "m_label = " << m_label << std::endl;

    for (const auto& lbl : labels) {
        //! Fields to be averaged
        amrex::Vector<std::string> fnames;
        std::string avg_type;
        const std::string pp_key = m_label + "." + lbl;
        std::cout << "lbl = " << lbl << std::endl;

        amrex::ParmParse pp1(pp_key);
        pp1.getarr("fields", fnames);
        ioutils::assert_with_message(
            ioutils::all_distinct(fnames),
            "Duplicates in " + pp_key + ".fields");
        pp1.get("averaging_type", avg_type);

        for (const auto& fname : fnames) {
            const std::string key = fname + "_" + avg_type;

            // Create the averaging entity
            m_averages.emplace_back(
                FieldTimeAverage::create(avg_type, m_sim, m_label, fname));

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
    if (found != m_registered.end()) {
        return m_registered[key]->average_field_name();
    }

    // Create and register new average
    m_averages.emplace_back(
        FieldTimeAverage::create(avg_type, m_sim, m_label, field_name));
    return m_averages.back()->average_field_name();
}

void TimeAveraging::post_advance_work()
{
    const auto& time = m_sim.time();
    const auto cur_time = time.new_time();
    const auto cur_step = time.time_index();
    const auto cur_dt = time.delta_t();

    // Check the following:
    //   1. if we are within the averaging time period requested by the user
    //   2. if we are phase averaging (i.e., only accumulating the average at a certain
    //      frequency), check to see if we are on a time step in which averaging should
    //      be performed.  This is useful for simulations with periodic behavior, such 
    //      as regular waves or actuator lines rotating at fixed rotor speed.
    const bool do_avg =
        (((cur_time >= m_start_time) && (cur_time < m_stop_time)) && 
         (cur_step % m_frequency == 0));
    if (!do_avg) {
        return;
    }

    const amrex::Real elapsed_time = (cur_time - m_start_time);
    m_last_avg_time = cur_time;
    std::cout << "Averaging " << m_label << std::endl;
    for (const auto& avg : m_averages) {
        std::cout << "Accumulating average.  Last accumulation done at: " << m_last_avg_time << ", elapsed time: " << elapsed_time << std::endl;
        (*avg)(time, m_filter, m_frequency, elapsed_time);
    }
}

} // namespace amr_wind::averaging
