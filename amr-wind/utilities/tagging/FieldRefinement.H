#ifndef FIELDREFINEMENT_H
#define FIELDREFINEMENT_H

#include "amr-wind/utilities/tagging/RefinementCriteria.H"

namespace amr_wind {
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
                            std::abs(farr(i + 1, j, k) - farr(i, j, k));
                        const amrex::Real ayp =
                            std::abs(farr(i, j + 1, k) - farr(i, j, k));
                        const amrex::Real azp =
                            std::abs(farr(i, j, k + 1) - farr(i, j, k));
                        const amrex::Real axm =
                            std::abs(farr(i - 1, j, k) - farr(i, j, k));
                        const amrex::Real aym =
                            std::abs(farr(i, j - 1, k) - farr(i, j, k));
                        const amrex::Real azm =
                            std::abs(farr(i, j, k - 1) - farr(i, j, k));
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

private:
    const CFDSim& m_sim;

    Field* m_field{nullptr};
    IntField* m_int_field{nullptr};

    amrex::Vector<amrex::Real> m_field_error;
    amrex::Vector<amrex::Real> m_grad_error;

    int m_max_lev_field{-1};
    int m_max_lev_grad{-1};
};

} // namespace amr_wind

#endif /* FIELDREFINEMENT_H */
