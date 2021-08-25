#include "amr-wind/utilities/tagging/TracerRefinement.H"
#include "amr-wind/CFDSim.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

TracerRefinement::TracerRefinement(const CFDSim& sim)
    : m_sim(sim)
    , m_max_lev(m_sim.mesh().maxLevel())
    , m_tracer_scalar(sim.repo().get_field("passive_tracer"))
{}

void TracerRefinement::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);
    pp.query("max_level", m_max_lev);
    pp.query("value", m_value);
}

void TracerRefinement::operator()(
    int level, amrex::TagBoxArray& tags, amrex::Real, int)
{
    if (level > m_max_lev) return;

    const amrex::Real value = m_value;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(m_tracer_scalar(level), amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const auto& bx = mfi.tilebox();
        const auto& tracer_arr = m_tracer_scalar(level).const_array(mfi);
        const auto& tag = tags.array(mfi);

        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                if (tracer_arr(i, j, k) >= value)
                    tag(i, j, k) = amrex::TagBox::SET;
            });
    }
}

} // namespace amr_wind
