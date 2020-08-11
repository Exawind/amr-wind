#include "amr-wind/utilities/PostProcessing.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {

PostProcessManager::PostProcessManager(CFDSim& sim) : m_sim(sim) {}

void PostProcessManager::initialize()
{
    amrex::Vector<std::string> pnames;
    amrex::ParmParse pp("incflo");
    pp.queryarr("post_processing", pnames);
    for (const auto& label: pnames) {
        std::string ptype = "Sampling";
        amrex::ParmParse pp1(label);
        pp1.query("type", ptype);
        m_post.emplace_back(PostProcessBase::create(ptype, m_sim, label));
    }

    for (auto& post: m_post) {
        post->initialize();
        post->post_advance_work();
    }
}

void PostProcessManager::post_advance_work()
{
    for (auto& post: m_post)
        post->post_advance_work();
}

} // namespace amr_wind
