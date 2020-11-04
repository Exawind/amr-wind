#include "amr-wind/utilities/tagging/QCriterionRefinement.H"
#include "amr-wind/CFDSim.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

QCriterionRefinement::QCriterionRefinement(const CFDSim& sim)
    : m_sim(sim)
    , m_qc_value(
          m_sim.mesh().maxLevel() + 1, std::numeric_limits<amrex::Real>::max())
{}

void QCriterionRefinement::initialize(const std::string& key)
{
    std::string fname = "velocity";

    const auto& repo = m_sim.repo();
    if (!repo.field_exists(fname)) {
        amrex::Abort("FieldRefinement: Cannot find field = " + fname);
    }
    m_vel = &(m_sim.repo().get_field(fname));

    amrex::Vector<double> qc_value;
    amrex::ParmParse pp(key);

    pp.queryarr("values", qc_value);

    if (qc_value.size() == 0u)
        amrex::Abort(
            "QCriterionRefinement: Must specify at least one of value");

    {
        size_t fcount = std::min(qc_value.size(), m_qc_value.size());
        for (size_t i = 0; i < fcount; ++i) m_qc_value[i] = qc_value[i];
        m_max_lev_field = fcount - 1;
    }
}

void QCriterionRefinement::operator()(
    int level, amrex::TagBoxArray& tags, amrex::Real time, int)
{
    const bool tag_field = level <= m_max_lev_field;

    m_vel->fillpatch(level, time, (*m_vel)(level), 1);

    const auto& mfab = (*m_vel)(level);
    const auto& idx = m_sim.repo().mesh().Geom(level).InvCellSizeArray();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mfab, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const auto& bx = mfi.tilebox();
        const auto& tag = tags.array(mfi);
        const auto& vel = mfab.const_array(mfi);

        if (tag_field) {
            const auto qc_val = m_qc_value[level];
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

            // this ignores wall stencils
#if 1
                    const auto ux =
                        0.5 * (vel(i + 1, j, k, 0) - vel(i - 1, j, k, 0)) *
                        idx[0];
                    const auto vx =
                        0.5 * (vel(i + 1, j, k, 1) - vel(i - 1, j, k, 1)) *
                        idx[0];
                    const auto wx =
                        0.5 * (vel(i + 1, j, k, 2) - vel(i - 1, j, k, 2)) *
                        idx[0];

                    const auto uy =
                        0.5 * (vel(i, j + 1, k, 0) - vel(i, j - 1, k, 0)) *
                        idx[1];
                    const auto vy =
                        0.5 * (vel(i, j + 1, k, 1) - vel(i, j - 1, k, 1)) *
                        idx[1];
                    const auto wy =
                        0.5 * (vel(i, j + 1, k, 2) - vel(i, j - 1, k, 2)) *
                        idx[1];

                    const auto uz =
                        0.5 * (vel(i, j, k + 1, 0) - vel(i, j, k - 1, 0)) *
                        idx[2];
                    const auto vz =
                        0.5 * (vel(i, j, k + 1, 1) - vel(i, j, k - 1, 1)) *
                        idx[2];
                    const auto wz =
                        0.5 * (vel(i, j, k + 1, 2) - vel(i, j, k - 1, 2)) *
                        idx[2];
#else
                    const auto ux = amrex::max(
                                        vel(i + 1, j, k, 0) - vel(i, j, k, 0),
                                        vel(i, j, k, 0) - vel(i - 1, j, k, 0)) *
                                    idx[0];
                    const auto vx = amrex::max(
                                        vel(i + 1, j, k, 1) - vel(i, j, k, 1),
                                        vel(i, j, k, 1) - vel(i - 1, j, k, 1)) *
                                    idx[0];
                    const auto wx = amrex::max(
                                        vel(i + 1, j, k, 2) - vel(i, j, k, 2),
                                        vel(i, j, k, 2) - vel(i - 1, j, k, 2)) *
                                    idx[0];

                    const auto uy = amrex::max(
                                        vel(i, j + 1, k, 0) - vel(i, j, k, 0),
                                        vel(i, j, k, 0) - vel(i, j - 1, k, 0)) *
                                    idx[1];
                    const auto vy = amrex::max(
                                        vel(i, j + 1, k, 1) - vel(i, j, k, 1),
                                        vel(i, j, k, 1) - vel(i, j - 1, k, 1)) *
                                    idx[1];
                    const auto wy = amrex::max(
                                        vel(i, j + 1, k, 2) - vel(i, j, k, 2),
                                        vel(i, j, k, 2) - vel(i, j - 1, k, 2)) *
                                    idx[1];

                    const auto uz = amrex::max(
                                        vel(i, j, k + 1, 0) - vel(i, j, k, 0),
                                        vel(i, j, k, 0) - vel(i, j, k - 1, 0)) *
                                    idx[2];
                    const auto vz = amrex::max(
                                        vel(i, j, k + 1, 1) - vel(i, j, k, 1),
                                        vel(i, j, k, 1) - vel(i, j, k - 1, 1)) *
                                    idx[2];
                    const auto wz = amrex::max(
                                        vel(i, j, k + 1, 2) - vel(i, j, k, 2),
                                        vel(i, j, k, 2) - vel(i, j, k - 1, 2)) *
                                    idx[2];
#endif
                    const auto S2 = 2.0 * ux * ux + 2.0 * vy * vy +
                                    2.0 * wz * wz + std::pow(uy + vx, 2) +
                                    std::pow(vz + wy, 2) + std::pow(wx + uz, 2);

                    const auto W2 = std::pow(uy - vx, 2) +
                                    std::pow(vz - wy, 2) + std::pow(wx - uz, 2);

                    const auto qc = 0.5 * (W2 - S2);

                    // absolute value is done on purpose and will tag a larger
                    // radius around vortices
                    if (amrex::Math::abs(qc) > qc_val)
                        tag(i, j, k) = amrex::TagBox::SET;
                });
        }
    }
}

} // namespace amr_wind
