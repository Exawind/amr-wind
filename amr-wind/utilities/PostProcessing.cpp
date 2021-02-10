#include "amr-wind/utilities/PostProcessing.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/averaging/TimeAveraging.H"

#include "AMReX_ParmParse.H"

#include <set>

namespace amr_wind {

namespace {
void perform_checks(
    std::set<std::string>& registered_types, const std::string ptype)
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

void PostProcessManager::initialize()
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
        post->initialize();
        post->post_advance_work();
    }
}

void PostProcessManager::post_advance_work()
{
    for (auto& post : m_post) post->post_advance_work();
}

} // namespace amr_wind
