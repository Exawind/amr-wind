#include "amr-wind/utilities/PostProcessing.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/averaging/TimeAveraging.H"

#include "AMReX_ParmParse.H"

#include <set>

namespace amr_wind {

namespace {
void perform_checks(
    std::set<std::string>& registered_types, const std::string& ptype)
{
    const auto found = registered_types.find(ptype);

    if (found == registered_types.end()) {
        registered_types.insert(ptype);
        return;
    }

    if (ptype == averaging::TimeAveraging::identifier()) {
        amrex::Abort(
            "PostProcessing: Multiple instances of TimeAveraging not allowed");
    }
}
} // namespace

PostProcessManager::PostProcessManager(CFDSim& sim) : m_sim(sim) {}

void PostProcessManager::pre_init_actions()
{
    amrex::Vector<std::string> pnames;
    amrex::ParmParse pp("incflo");
    pp.queryarr("post_processing", pnames);
    std::set<std::string> registered_types;

    for (const auto& label : pnames) {
        std::string ptype = "Sampling";
        amrex::ParmParse pp1(label);
        pp1.query("type", ptype);

        perform_checks(registered_types, ptype);
        m_post.emplace_back(PostProcessBase::create(ptype, m_sim, label));
    }

    for (auto& post : m_post) {
        post->pre_init_actions();
    }
}

void PostProcessManager::post_init_actions()
{
    for (auto& post : m_post) {
        post->initialize();
        m_sim.time().add_postproc_dt_parameters(
            post->enforce_dt(), post->enforce_dt_tolerance(),
            post->output_time_interval(), post->output_time_delay());
        post->post_advance_work();
    }
    // Calculate and get minimum tolerance
    m_sim.time().calculate_minimum_enforce_dt_abs_tol();
    auto tol = m_sim.time().get_minimum_enforce_dt_abs_tol();
    for (auto& post : m_post) {
        if (post->do_output_now(
                m_sim.time().time_index(), m_sim.time().new_time(),
                m_sim.time().delta_t(), tol)) {
            post->output_actions();
        }
    }
}

void PostProcessManager::post_advance_work()
{
    // Calculate and get minimum tolerance
    m_sim.time().calculate_minimum_enforce_dt_abs_tol();
    auto tol = m_sim.time().get_minimum_enforce_dt_abs_tol();
    for (auto& post : m_post) {
        post->post_advance_work();
        if (post->do_output_now(
                m_sim.time().time_index(), m_sim.time().new_time(),
                m_sim.time().delta_t(), tol)) {
            post->output_actions();
        }
    }
}

void PostProcessManager::post_regrid_actions()
{
    for (auto& post : m_post) {
        post->post_regrid_actions();
    }
}

} // namespace amr_wind
