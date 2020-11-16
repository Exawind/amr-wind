#include "amr-wind/utilities/tagging/OversetRefinement.H"
#include "amr-wind/CFDSim.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

OversetRefinement::OversetRefinement(const CFDSim& sim)
    : m_sim(sim), m_max_lev(m_sim.mesh().maxLevel())
{}

void OversetRefinement::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);
    pp.query("max_level", m_max_lev);
}

void OversetRefinement::operator()(
    int level, amrex::TagBoxArray& tags, amrex::Real, int)
{
    if (level > m_max_lev) return;

    const auto& ibcell = m_sim.repo().get_int_field("iblank_cell");
    const auto& ibfab = ibcell(level);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(ibfab, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const auto& bx = mfi.tilebox();
        const auto& ibarr = ibfab.const_array(mfi);
        const auto& tag = tags.array(mfi);

        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const int axp =
                    amrex::Math::abs(ibarr(i + 1, j, k) - ibarr(i, j, k));
                const int ayp =
                    amrex::Math::abs(ibarr(i, j + 1, k) - ibarr(i, j, k));
                const int azp =
                    amrex::Math::abs(ibarr(i, j, k + 1) - ibarr(i, j, k));
                const int axm =
                    amrex::Math::abs(ibarr(i - 1, j, k) - ibarr(i, j, k));
                const int aym =
                    amrex::Math::abs(ibarr(i, j - 1, k) - ibarr(i, j, k));
                const int azm =
                    amrex::Math::abs(ibarr(i, j, k - 1) - ibarr(i, j, k));
                const int ax = amrex::max(axp, axm);
                const int ay = amrex::max(ayp, aym);
                const int az = amrex::max(azp, azm);
                if (amrex::max(ax, ay, az) > 1)
                    tag(i, j, k) = amrex::TagBox::SET;
            });
    }
}

} // namespace amr_wind
