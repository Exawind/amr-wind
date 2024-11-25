#include "amr-wind/physics/ForestDrag.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParReduce.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/utilities/index_operations.H"

namespace amr_wind::forestdrag {

namespace {} // namespace

ForestDrag::ForestDrag(CFDSim& sim)
    : m_sim(sim)
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_forest_drag(sim.repo().declare_field("forest_drag", 1, 1, 1))
    , m_forest_blank(sim.repo().declare_field("forest_blank", 1, 1, 1))
{
    std::string forestfile("forest.amrwind");
    amrex::ParmParse pp(identifier());
    pp.query("forest_file", forestfile);
    std::ifstream file(forestfile, std::ios::in);
    if (!file.good()) {
        amrex::Abort("Cannot find file " + forestfile);
    }
    //! TreeType xc yc height diameter cd lai laimax
    amrex::Real value1, value2, value3, value4, value5, value6, value7, value8;
    while (file >> value1 >> value2 >> value3 >> value4 >> value5 >> value6 >>
           value7 >> value8) {
        Forest f;
        f.m_type_forest = value1;
        f.m_x_forest = value2;
        f.m_y_forest = value3;
        f.m_height_forest = value4;
        f.m_diameter_forest = value5;
        f.m_cd_forest = value6;
        f.m_lai_forest = value7;
        f.m_laimax_forest = value8;
        m_forests.push_back(f);
    }
    file.close();
    m_d_forests.resize(m_forests.size());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_forests.begin(), m_forests.end(),
        m_d_forests.begin());
    m_sim.io_manager().register_output_var("forest_drag");
    m_sim.io_manager().register_output_var("forest_blank");
}

void ForestDrag::post_init_actions()
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::post_init_actions");
    const auto& geom_vec = m_sim.repo().mesh().Geom();
    const int nlevels = m_sim.repo().num_active_levels();
    for (int level = 0; level < nlevels; ++level) {
        const auto& geom = geom_vec[level];
        const auto& dx = geom.CellSizeArray();
        const auto& prob_lo = geom.ProbLoArray();
        auto& velocity = m_velocity(level);
        auto& drag = m_forest_drag(level);
        auto& blank = m_forest_blank(level);
        const auto* forests_ptr = m_d_forests.data();
        const auto forestSize = m_d_forests.size();
        for (unsigned ii = 0; ii < forestSize; ++ii) {
            const amrex::Real treelaimax = forests_ptr[ii].calc_lm();
            const amrex::Real x_forest = forests_ptr[ii].m_x_forest;
            const amrex::Real y_forest = forests_ptr[ii].m_y_forest;
            const amrex::Real height_forest = forests_ptr[ii].m_height_forest;
            const amrex::Real diameter_forest =
                forests_ptr[ii].m_diameter_forest;
            const amrex::Real cd_forest = forests_ptr[ii].m_cd_forest;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
                const auto& vbx = mfi.validbox();
                auto levelDrag = drag.array(mfi);
                auto levelBlank = blank.array(mfi);
                const amrex::Real x1 = x_forest - 0.5 * diameter_forest;
                const amrex::Real y1 = y_forest - 0.5 * diameter_forest;
                const amrex::Real z1 = prob_lo[2];
                const amrex::Real x2 = x_forest + 0.5 * diameter_forest;
                const amrex::Real y2 = y_forest + 0.5 * diameter_forest;
                const amrex::Real z2 = height_forest;
                const amrex::RealBox& real_treebox{x1, y1, z1, x2, y2, z2};
                const auto valid_treebox =
                    amr_wind::utils::realbox_to_box(real_treebox, geom);
                if (vbx.contains(valid_treebox)) {
                    amrex::ParallelFor(
                        vbx,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            // compute the source term
                            const amrex::Real x =
                                prob_lo[0] + (i + 0.5) * dx[0];
                            const amrex::Real y =
                                prob_lo[1] + (j + 0.5) * dx[1];
                            const amrex::Real z =
                                prob_lo[2] + (k + 0.5) * dx[2];
                            const amrex::Real radius = std::sqrt(
                                (x - x_forest) * (x - x_forest) +
                                (y - y_forest) * (y - y_forest));
                            if (z <= height_forest &&
                                radius <= (0.5 * diameter_forest)) {
                                levelBlank(i, j, k) = ii + 1;
                                levelDrag(i, j, k) =
                                    cd_forest *
                                    forests_ptr[ii].calc_af(z, treelaimax);
                            }
                        });
                }
            }
        }
    }
}

void ForestDrag::pre_init_actions()
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::pre_init_actions");
}
} // namespace amr_wind::forestdrag
