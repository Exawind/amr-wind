#ifndef FIELDREFINEMENT_H
#define FIELDREFINEMENT_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tagging/RefinementCriteria.H"

namespace amr_wind {
class CFDSim;
class Field;
class IntField;

/** AMR refinement using a given field (e.g., density)
 *  \ingroup amr_utils
 *
 *  ```
 *  tagging.labels = t1
 *  tagging/t1.type = FieldRefinement
 *  tagging/t1.field_name = density
 *  # Tagging based on field value
 *  # tagging/t1.field_error = 10.0 10.0

 *  # Tagging based on local gradient
 *  tagging/t1.grad_error = 0.1 0.1 0.1
 *  ```
 */
class FieldRefinement : public RefinementCriteria::Register<FieldRefinement>
{
public:
    static std::string identifier() { return "FieldRefinement"; }

    explicit FieldRefinement(const CFDSim& sim);

    ~FieldRefinement() override = default;

    //! Read input file and initialize boxarray used to refine each level
    void initialize(const std::string& key) override;

    void operator()(
        const int level,
        amrex::TagBoxArray& tags,
        const amrex::Real time,
        const int ngrow) override;

    template <typename MF>
    void tag(const int level, amrex::TagBoxArray& tags, const MF& mfab)
    {
        const bool tag_field = level <= m_max_lev_field;
        const bool tag_grad = level <= m_max_lev_grad;
        const auto& geom = m_sim.repo().mesh().Geom(level);
        const auto& prob_lo = geom.ProbLoArray();
        const auto& dx = geom.CellSizeArray();
        const auto tagging_box = m_tagging_box;

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
                        const amrex::RealVect coord = {AMREX_D_DECL(
                            prob_lo[0] + (i + 0.5) * dx[0],
                            prob_lo[1] + (j + 0.5) * dx[1],
                            prob_lo[2] + (k + 0.5) * dx[2])};

                        if ((farr(i, j, k) > fld_err) &&
                            (tagging_box.contains(coord))) {
                            tag(i, j, k) = amrex::TagBox::SET;
                        }
                    });
            }

            if (tag_grad) {
                const auto gerr = m_grad_error[level];
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::RealVect coord = {AMREX_D_DECL(
                            prob_lo[0] + (i + 0.5) * dx[0],
                            prob_lo[1] + (j + 0.5) * dx[1],
                            prob_lo[2] + (k + 0.5) * dx[2])};
                        const auto axp =
                            std::abs(farr(i + 1, j, k) - farr(i, j, k));
                        const auto ayp =
                            std::abs(farr(i, j + 1, k) - farr(i, j, k));
                        const auto azp =
                            std::abs(farr(i, j, k + 1) - farr(i, j, k));
                        const auto axm =
                            std::abs(farr(i - 1, j, k) - farr(i, j, k));
                        const auto aym =
                            std::abs(farr(i, j - 1, k) - farr(i, j, k));
                        const auto azm =
                            std::abs(farr(i, j, k - 1) - farr(i, j, k));
                        const auto ax = amrex::max(axp, axm);
                        const auto ay = amrex::max(ayp, aym);
                        const auto az = amrex::max(azp, azm);
                        if ((amrex::max(ax, ay, az) >= gerr) &&
                            (tagging_box.contains(coord))) {
                            tag(i, j, k) = amrex::TagBox::SET;
                        }
                    });
            }
        }
    }

private:
    const CFDSim& m_sim;

    Field* m_field{nullptr};
    IntField* m_int_field{nullptr};

    amrex::Vector<amrex::Real> m_field_error;
    amrex::Vector<amrex::Real> m_grad_error;

    int m_max_lev_field{-1};
    int m_max_lev_grad{-1};
    amrex::RealBox m_tagging_box;
};

} // namespace amr_wind

#endif /* FIELDREFINEMENT_H */
