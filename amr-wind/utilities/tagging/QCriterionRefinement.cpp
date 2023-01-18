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

    if (qc_value.empty()) {
        amrex::Abort(
            "QCriterionRefinement: Must specify at least one of value");
    }

    {
        const int fcount = std::min(
            static_cast<int>(qc_value.size()),
            static_cast<int>(m_qc_value.size()));
        for (int i = 0; i < fcount; ++i) {
            m_qc_value[i] = qc_value[i];
        }
        m_max_lev_field = fcount - 1;
    }

    pp.query("nondim", m_nondim);
}

void QCriterionRefinement::operator()(
    int level, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/)
{
    const bool tag_field = level <= m_max_lev_field;

    if (!tag_field) {
        return;
    }

    m_vel->fillpatch(level, time, (*m_vel)(level), 1);

    const auto& mfab = (*m_vel)(level);
    const auto& idx = m_sim.repo().mesh().Geom(level).InvCellSizeArray();

    const auto nondim = m_nondim;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mfab, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const auto& bx = mfi.tilebox();
        const auto& tag = tags.array(mfi);
        const auto& vel = mfab.const_array(mfi);
        const auto qc_val = m_qc_value[level];

        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // TODO: ignoring wall stencils for now
                const auto ux =
                    0.5 * (vel(i + 1, j, k, 0) - vel(i - 1, j, k, 0)) * idx[0];
                const auto vx =
                    0.5 * (vel(i + 1, j, k, 1) - vel(i - 1, j, k, 1)) * idx[0];
                const auto wx =
                    0.5 * (vel(i + 1, j, k, 2) - vel(i - 1, j, k, 2)) * idx[0];

                const auto uy =
                    0.5 * (vel(i, j + 1, k, 0) - vel(i, j - 1, k, 0)) * idx[1];
                const auto vy =
                    0.5 * (vel(i, j + 1, k, 1) - vel(i, j - 1, k, 1)) * idx[1];
                const auto wy =
                    0.5 * (vel(i, j + 1, k, 2) - vel(i, j - 1, k, 2)) * idx[1];

                const auto uz =
                    0.5 * (vel(i, j, k + 1, 0) - vel(i, j, k - 1, 0)) * idx[2];
                const auto vz =
                    0.5 * (vel(i, j, k + 1, 1) - vel(i, j, k - 1, 1)) * idx[2];
                const auto wz =
                    0.5 * (vel(i, j, k + 1, 2) - vel(i, j, k - 1, 2)) * idx[2];

                const auto S2 =
                    ux * ux + vy * vy + wz * wz + 0.5 * std::pow(uy + vx, 2) +
                    0.5 * std::pow(vz + wy, 2) + 0.5 * std::pow(wx + uz, 2);

                const auto W2 = 0.5 * std::pow(uy - vx, 2) +
                                0.5 * std::pow(vz - wy, 2) +
                                0.5 * std::pow(wx - uz, 2);

                const auto qc = 0.5 * (W2 - S2);
                const auto qc_nondim =
                    0.5 * (W2 / amrex::max(S2, 1.0e-12) - 1.0);

                if ((nondim && qc_nondim > qc_val) ||
                    (!nondim && amrex::Math::abs(qc) > qc_val)) {
                    tag(i, j, k) = amrex::TagBox::SET;
                }
            });
    }
}

} // namespace amr_wind
