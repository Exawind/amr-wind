#include "amr-wind/physics/ForestDrag.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParReduce.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/IOManager.H"

namespace amr_wind::forestdrag {

ForestDrag::ForestDrag(CFDSim& sim)
    : m_sim(sim)
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

        // Only keep a forest if the ranks owns it
        const int nlevels = m_sim.repo().num_active_levels();
        for (int level = 0; level < nlevels; ++level) {
            const auto& geom = m_sim.repo().mesh().Geom(level);
            const auto bx = f.bounding_box(geom);
            const auto& ba = m_sim.repo().mesh().boxArray(level);
            if (ba.contains(bx)) {
                m_forests.push_back(f);
                break;
            }
        }
    }
    file.close();
    m_d_forests.resize(m_forests.size());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_forests.begin(), m_forests.end(),
        m_d_forests.begin());
    m_sim.io_manager().register_output_var("forest_drag");
    m_sim.io_manager().register_output_var("forest_blank");

    m_forest_drag.setVal(0.0);
    m_forest_blank.setVal(-1.0);
    m_forest_blank.set_default_fillpatch_bc(m_sim.time());
    m_forest_drag.set_default_fillpatch_bc(m_sim.time());
}

void ForestDrag::initialize_fields(int level, const amrex::Geometry& geom)
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::initialize_fields");
    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    auto& drag = m_forest_drag(level);
    auto& blank = m_forest_blank(level);
    drag.setVal(0.0);
    blank.setVal(-1.0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(m_forest_drag(level)); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.tilebox();
        for (int nf = 0; nf < static_cast<int>(m_forests.size()); nf++) {
            if (vbx.contains(m_forests[nf].bounding_box(geom))) {
                const auto& levelDrag = drag.array(mfi);
                const auto& levelBlank = blank.array(mfi);
                const auto* d_forests = m_d_forests.data();
                amrex::ParallelFor(
                    vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const auto x = prob_lo[0] + (i + 0.5) * dx[0];
                        const auto y = prob_lo[1] + (j + 0.5) * dx[1];
                        const auto z = prob_lo[2] + (k + 0.5) * dx[2];
                        const auto& fst = d_forests[nf];
                        const auto radius = std::sqrt(
                            (x - fst.m_x_forest) * (x - fst.m_x_forest) +
                            (y - fst.m_y_forest) * (y - fst.m_y_forest));
                        if (z <= fst.m_height_forest &&
                            radius <= (0.5 * fst.m_diameter_forest)) {
                            const auto treelaimax = fst.calc_lm();
                            levelBlank(i, j, k) = nf;
                            levelDrag(i, j, k) +=
                                fst.m_cd_forest * fst.calc_af(z, treelaimax);
                        }
                    });
            }
        }
    }
}

void ForestDrag::post_regrid_actions()
{
    const int nlevels = m_sim.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        initialize_fields(lev, m_sim.repo().mesh().Geom(lev));
    }
}
} // namespace amr_wind::forestdrag
