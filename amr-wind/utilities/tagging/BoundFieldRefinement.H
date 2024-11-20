#ifndef BOUNDFIELDREFINEMENT_H
#define BOUNDFIELDREFINEMENT_H

#include "amr-wind/utilities/tagging/RefinementCriteria.H"

namespace amr_wind {
class Field;
class IntField;

/** AMR refinement using a given field (e.g., density)
 *  \ingroup amr_utils
 *
 *  ```
 *  tagging.labels = t1
 *  tagging/t1.type = BoundFieldRefinement
 *  tagging/t1.field_name = density
 *  # Tagging based on field value
 *  # tagging/t1.field_error = 10.0 10.0

 *  # Tagging based on local gradient
 *  tagging/t1.grad_error = 0.1 0.1 0.1
 *  ```
 * Adds an option to define bounds for different level of refinement
 * from a user specified file
 * A sample file looks like this
 * xmin xmax ymin ymax level
 * -500 500 -500 500 1
 * -200 200 -200 200 2
 * This helps to provide higher refinement levels near steep cliffs
 * and keep refinement levels lower away from region of interest
 */
class BoundFieldRefinement
    : public RefinementCriteria::Register<BoundFieldRefinement>
{
public:
    static std::string identifier() { return "BoundFieldRefinement"; }

    explicit BoundFieldRefinement(const CFDSim& sim);

    ~BoundFieldRefinement() override = default;

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
        const auto& geom = m_mesh.Geom(level);
        const auto& dx = geom.CellSizeArray();
        const auto& prob_lo = geom.ProbLoArray();
        const unsigned vsize = m_xmin_list.size();
        amrex::Gpu::DeviceVector<amrex::Real> m_device_xmin_list;
        amrex::Gpu::DeviceVector<amrex::Real> m_device_xmax_list;
        amrex::Gpu::DeviceVector<amrex::Real> m_device_ymin_list;
        amrex::Gpu::DeviceVector<amrex::Real> m_device_ymax_list;
        amrex::Gpu::DeviceVector<amrex::Real> m_device_level_list;
        m_device_xmin_list.resize(vsize);
        m_device_xmax_list.resize(vsize);
        m_device_ymin_list.resize(vsize);
        m_device_ymax_list.resize(vsize);
        m_device_level_list.resize(vsize);
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_xmin_list.begin(), m_xmin_list.end(),
            m_device_xmin_list.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_xmax_list.begin(), m_xmax_list.end(),
            m_device_xmax_list.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_ymin_list.begin(), m_ymin_list.end(),
            m_device_ymin_list.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_ymax_list.begin(), m_ymax_list.end(),
            m_device_ymax_list.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_level_list.begin(), m_level_list.end(),
            m_device_level_list.begin());
        const auto* device_xmin_list = m_device_xmin_list.data();
        const auto* device_xmax_list = m_device_xmax_list.data();
        const auto* device_ymin_list = m_device_ymin_list.data();
        const auto* device_ymax_list = m_device_ymax_list.data();
        const auto* device_level_list = m_device_level_list.data();
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
                        const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                        const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                        const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
                        bool add_cell = false;
                        for (unsigned ii = 0; ii < vsize; ++ii) {
                            if (x >= device_xmin_list[i] &&
                                x <= device_xmax_list[i] &&
                                y >= device_ymin_list[i] &&
                                y <= device_ymax_list[i] &&
                                level <= device_level_list[i] && z > dx[2]) {
                                add_cell = true;
                                break;
                            }
                        }
                        if (farr(i, j, k) > fld_err && add_cell) {
                            tag(i, j, k) = amrex::TagBox::SET;
                        }
                    });
            }
            if (tag_grad) {
                const auto gerr = m_grad_error[level];
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
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
                        const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                        const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                        const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
                        bool add_cell = false;
                        for (unsigned ii = 0; ii < vsize; ++ii) {
                            if (x >= device_xmin_list[i] &&
                                x <= device_xmax_list[i] &&
                                y >= device_ymin_list[i] &&
                                y <= device_ymax_list[i] &&
                                level <= device_level_list[i] && z > dx[2]) {
                                add_cell = true;
                                break;
                            }
                        }
                        if (amrex::max(ax, ay, az) >= gerr && add_cell) {
                            tag(i, j, k) = amrex::TagBox::SET;
                        }
                    });
            }
        }
    }

private:
    const CFDSim& m_sim;
    const amrex::AmrCore& m_mesh;
    Field* m_field{nullptr};
    IntField* m_int_field{nullptr};

    amrex::Vector<amrex::Real> m_field_error;
    amrex::Vector<amrex::Real> m_grad_error;
    std::string m_boundsfile;
    amrex::Vector<amrex::Real> m_xmin_list;
    amrex::Vector<amrex::Real> m_xmax_list;
    amrex::Vector<amrex::Real> m_ymin_list;
    amrex::Vector<amrex::Real> m_ymax_list;
    amrex::Vector<amrex::Real> m_level_list;
    int m_max_lev_field{-1};
    int m_max_lev_grad{-1};
};

} // namespace amr_wind

#endif /* BOUNDFIELDREFINEMENT_H */
