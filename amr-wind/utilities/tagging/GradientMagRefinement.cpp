#include "amr-wind/utilities/tagging/GradientMagRefinement.H"
#include "amr-wind/CFDSim.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

GradientMagRefinement::GradientMagRefinement(const CFDSim& sim)
    : m_sim(sim)
    , m_gradmag_value(
          m_sim.mesh().maxLevel() + 1, std::numeric_limits<amrex::Real>::max())
{}

void GradientMagRefinement::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);
    std::string fname;
    pp.query("field_name", fname);

    const auto& repo = m_sim.repo();
    if (!repo.field_exists(fname)) {
        amrex::Abort("GradientMagRefinement: Cannot find field = " + fname);
    }
    m_field = &(m_sim.repo().get_field(fname));

    amrex::Vector<double> gradmag_value;
    pp.queryarr("values", gradmag_value);

    if ((gradmag_value.empty())) {
        amrex::Abort("GradientMagRefinement: Must specify at least one value");
    }

    {
        const int fcount = std::min(
            static_cast<int>(gradmag_value.size()),
            static_cast<int>(m_gradmag_value.size()));
        for (int i = 0; i < fcount; ++i) {
            m_gradmag_value[i] = gradmag_value[i];
        }
        m_max_lev_field = fcount - 1;
    }
}

void GradientMagRefinement::operator()(
    int level, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/)
{
    const bool tag_field = level <= m_max_lev_field;
    if (tag_field) {
        m_field->fillpatch(level, time, (*m_field)(level), 1);
    }

    const auto& mfab = (*m_field)(level);
    const auto& idx = m_sim.repo().mesh().Geom(level).InvCellSizeArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mfab, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const auto& bx = mfi.tilebox();
        const auto& tag = tags.array(mfi);
        const auto& farr = mfab.const_array(mfi);
        const auto gradmag_val = m_gradmag_value[level];

        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // TODO: ignoring wall stencils for now

                const auto gx =
                    0.5 * (farr(i + 1, j, k) - farr(i - 1, j, k)) * idx[0];
                const auto gy =
                    0.5 * (farr(i, j + 1, k) - farr(i, j - 1, k)) * idx[1];
                const auto gz =
                    0.5 * (farr(i, j, k + 1) - farr(i, j, k - 1)) * idx[2];

                const auto grad_mag = sqrt(gx * gx + gy * gy + gz * gz);
                if (grad_mag > gradmag_val) {
                    tag(i, j, k) = amrex::TagBox::SET;
                }
            });
    }
}

} // namespace amr_wind
