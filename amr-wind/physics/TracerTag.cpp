#include "amr-wind/CFDSim.H"
#include "amr-wind/physics/TracerTag.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

TracerTag::TracerTag(CFDSim& sim) : m_sim(sim)
{
    auto& pseqn = sim.pde_manager().register_transport_pde("PassiveScalar");
    m_tracer = &(pseqn.fields().field);

    amrex::ParmParse pp("TracerTag");
    pp.query("act_src_threshold", m_src_threshold);
}

void TracerTag::initialize_fields(int level, const amrex::Geometry&)
{
    (*m_tracer)(level).setVal(0.0);
}

void TracerTag::post_init_actions()
{

    if (m_sim.repo().field_exists("actuator_src_term"))
        m_act_src = &(m_sim.repo().get_field("actuator_src_term"));

    if (m_sim.repo().field_exists("iblank_cell"))
        m_iblank = &(m_sim.repo().get_int_field("iblank_cell"));
}

void TracerTag::post_advance_work()
{

    if (!m_act_src && !m_iblank) {
        amrex::Print() << "Warning TracerTag activated but neither actuators "
                          "or overset are being used"
                       << std::endl;
        return;
    }

    for (int level = 0; level <= m_sim.mesh().finestLevel(); ++level) {
        for (amrex::MFIter mfi((*m_tracer)(level)); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.tilebox();
            auto tracer = (*m_tracer)(level).array(mfi);

            if (m_act_src) {

                auto src = (*m_act_src)(level).array(mfi);
                const amrex::Real src_threshold = m_src_threshold;

                amrex::ParallelFor(
                    vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real srcmag = std::sqrt(
                            src(i, j, k, 0) * src(i, j, k, 0) +
                            src(i, j, k, 1) * src(i, j, k, 1) +
                            src(i, j, k, 2) * src(i, j, k, 2));

                        if (srcmag > src_threshold) tracer(i, j, k) = 1.0;
                    });
            }

            if (m_iblank) {
                auto iblank = (*m_iblank)(level).array(mfi);
                const bool tag_fringe = m_tag_fringe;
                const bool tag_hole = m_tag_hole;
                amrex::ParallelFor(
                    vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        if (tag_fringe && iblank(i, j, k) == -1)
                            tracer(i, j, k) = 1.0;
                        if (tag_hole && iblank(i, j, k) == 0)
                            tracer(i, j, k) = 1.0;
                    });
            }
        }
    }
}

} // namespace amr_wind
