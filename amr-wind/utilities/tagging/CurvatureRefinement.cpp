#include "amr-wind/utilities/tagging/CurvatureRefinement.H"
#include "amr-wind/CFDSim.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

CurvatureRefinement::CurvatureRefinement(const CFDSim& sim)
    : m_sim(sim)
    , m_curv_value(
          m_sim.mesh().maxLevel() + 1, std::numeric_limits<amrex::Real>::max())
{}

void CurvatureRefinement::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);
    std::string fname;
    pp.query("field_name", fname);

    const auto& repo = m_sim.repo();
    if (!repo.field_exists(fname)) {
        amrex::Abort("CurvatureRefinement: Cannot find field = " + fname);
    }
    m_field = &(m_sim.repo().get_field(fname));

    amrex::Vector<double> curv_value;
    pp.queryarr("values", curv_value);

    if ((curv_value.empty())) {
        amrex::Abort("CurvatureRefinement: Must specify at least one value");
    }

    {
        const int fcount = std::min(
            static_cast<int>(curv_value.size()),
            static_cast<int>(m_curv_value.size()));
        for (int i = 0; i < fcount; ++i) {
            m_curv_value[i] = curv_value[i];
        }
        m_max_lev_field = fcount - 1;
    }
}

void CurvatureRefinement::operator()(
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
        const auto curv_val = m_curv_value[level];

        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // TODO: ignoring wall stencils for now

                const auto phixx = (farr(i + 1, j, k) - 2.0 * farr(i, j, k) +
                                    farr(i - 1, j, k)) *
                                   idx[0] * idx[0];
                const auto phiyy = (farr(i, j + 1, k) - 2.0 * farr(i, j, k) +
                                    farr(i, j - 1, k)) *
                                   idx[0] * idx[0];
                const auto phizz = (farr(i, j, k + 1) - 2.0 * farr(i, j, k) +
                                    farr(i, j, k - 1)) *
                                   idx[0] * idx[0];

                const auto phiz =
                    0.5 * (farr(i, j, k + 1) - farr(i, j, k - 1)) * idx[2];
                const auto phiz_ip1 =
                    0.5 * (farr(i + 1, j, k + 1) - farr(i + 1, j, k - 1)) *
                    idx[2];
                const auto phiz_im1 =
                    0.5 * (farr(i - 1, j, k + 1) - farr(i - 1, j, k - 1)) *
                    idx[2];
                const auto phiz_jp1 =
                    0.5 * (farr(i, j + 1, k + 1) - farr(i, j + 1, k - 1)) *
                    idx[2];
                const auto phiz_jm1 =
                    0.5 * (farr(i, j - 1, k + 1) - farr(i, j - 1, k - 1)) *
                    idx[2];

                const auto phiy =
                    0.5 * (farr(i, j + 1, k) - farr(i, j - 1, k)) * idx[1];
                const auto phiy_ip1 =
                    0.5 * (farr(i + 1, j + 1, k) - farr(i + 1, j - 1, k)) *
                    idx[1];
                const auto phiy_im1 =
                    0.5 * (farr(i - 1, j + 1, k) - farr(i - 1, j - 1, k)) *
                    idx[1];
                const auto phiyz = 0.5 * (phiz_jp1 - phiz_jm1) * idx[1];

                const auto phix =
                    0.5 * (farr(i + 1, j, k) - farr(i - 1, j, k)) * idx[0];
                const auto phixy = 0.5 * (phiy_ip1 - phiy_im1) * idx[0];
                const auto phixz = 0.5 * (phiz_ip1 - phiz_im1) * idx[0];

                const auto curv_mag =
                    std::abs(
                        phix * phix * phiyy - 2. * phix * phiy * phixy +
                        phiy * phiy * phixx + phix * phix * phizz -
                        2. * phix * phiz * phixz + phiz * phiz * phixx +
                        phiy * phiy * phizz - 2. * phiy * phiz * phiyz +
                        phiz * phiz * phiyy) /
                    std::pow(phix * phix + phiy * phiy + phiz * phiz, 1.5);
                const auto curv_min =
                    std::min(curv_val, std::cbrt(idx[0] * idx[1] * idx[2]));
                if (curv_mag > curv_min) {
                    tag(i, j, k) = amrex::TagBox::SET;
                }
            });
    }
}

} // namespace amr_wind
