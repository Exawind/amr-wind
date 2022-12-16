#include "amr-wind/utilities/tagging/FieldRefinement.H"
#include "amr-wind/CFDSim.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

FieldRefinement::FieldRefinement(const CFDSim& sim)
    : m_sim(sim)
    , m_field_error(
          m_sim.mesh().maxLevel() + 1, std::numeric_limits<amrex::Real>::max())
    , m_grad_error(
          m_sim.mesh().maxLevel() + 1, std::numeric_limits<amrex::Real>::max())
{}

void FieldRefinement::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);
    std::string fname;
    pp.query("field_name", fname);

    const auto& repo = m_sim.repo();
    if (!repo.field_exists(fname)) {
        amrex::Abort("FieldRefinement: Cannot find field = " + fname);
    }
    m_field = &(m_sim.repo().get_field(fname));

    amrex::Vector<double> field_err;
    amrex::Vector<double> grad_err;
    pp.queryarr("field_error", field_err);
    pp.queryarr("grad_error", grad_err);

    if ((field_err.empty()) && (grad_err.empty())) {
        amrex::Abort(
            "FieldRefinement: Must specify at least one of field_error or "
            "grad_error");
    }

    {
        const int fcount = std::min(
            static_cast<int>(field_err.size()),
            static_cast<int>(m_field_error.size()));
        for (int i = 0; i < fcount; ++i) {
            m_field_error[i] = field_err[i];
        }
        m_max_lev_field = fcount - 1;
    }
    {
        const int fcount = std::min(
            static_cast<int>(grad_err.size()),
            static_cast<int>(m_grad_error.size()));
        for (int i = 0; i < fcount; ++i) {
            m_grad_error[i] = grad_err[i];
        }
        m_max_lev_grad = fcount - 1;
    }
}

void FieldRefinement::operator()(
    int level, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/)
{
    const bool tag_field = level <= m_max_lev_field;
    const bool tag_grad = level <= m_max_lev_grad;
    if (tag_grad) {
        m_field->fillpatch(level, time, (*m_field)(level), 1);
    }

    const auto& mfab = (*m_field)(level);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mfab, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const auto& bx = mfi.tilebox();
        const auto& tag = tags.array(mfi);
        const auto& farr = mfab.const_array(mfi);

        if (tag_field) {
            const auto fld_err = m_field_error[level];
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    if (farr(i, j, k) > fld_err) {
                        tag(i, j, k) = amrex::TagBox::SET;
                    }
                });
        }

        if (tag_grad) {
            const auto gerr = m_grad_error[level];
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real axp =
                        amrex::Math::abs(farr(i + 1, j, k) - farr(i, j, k));
                    const amrex::Real ayp =
                        amrex::Math::abs(farr(i, j + 1, k) - farr(i, j, k));
                    const amrex::Real azp =
                        amrex::Math::abs(farr(i, j, k + 1) - farr(i, j, k));
                    const amrex::Real axm =
                        amrex::Math::abs(farr(i - 1, j, k) - farr(i, j, k));
                    const amrex::Real aym =
                        amrex::Math::abs(farr(i, j - 1, k) - farr(i, j, k));
                    const amrex::Real azm =
                        amrex::Math::abs(farr(i, j, k - 1) - farr(i, j, k));
                    const amrex::Real ax = amrex::max(axp, axm);
                    const amrex::Real ay = amrex::max(ayp, aym);
                    const amrex::Real az = amrex::max(azp, azm);
                    if (amrex::max(ax, ay, az) >= gerr) {
                        tag(i, j, k) = amrex::TagBox::SET;
                    }
                });
        }
    }
}

} // namespace amr_wind
