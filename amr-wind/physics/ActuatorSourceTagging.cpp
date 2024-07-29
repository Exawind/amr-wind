#include "amr-wind/CFDSim.H"
#include "amr-wind/physics/ActuatorSourceTagging.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

ActuatorSourceTagging::ActuatorSourceTagging(CFDSim& sim) : m_sim(sim)
{
    auto& pseqn = sim.pde_manager().register_transport_pde("PassiveScalar");
    m_tracer = &(pseqn.fields().field);

    amrex::ParmParse pp("ActuatorSourceTagging");
    pp.query("act_src_threshold", m_src_threshold);
}

void ActuatorSourceTagging::initialize_fields(int level, const amrex::Geometry&)
{
    (*m_tracer)(level).setVal(0.0);
}

void ActuatorSourceTagging::post_init_actions()
{

    if (m_sim.repo().field_exists("actuator_src_term"))
        m_act_src = &(m_sim.repo().get_field("actuator_src_term"));

    if (m_sim.repo().field_exists("iblank_cell"))
        m_iblank = &(m_sim.repo().get_int_field("iblank_cell"));
}

void ActuatorSourceTagging::post_advance_work()
{

    if (!m_act_src && !m_iblank) {
        amrex::Print()
            << "Warning ActuatorSourceTagging activated but neither actuators "
               "or overset are being used"
            << std::endl;
        return;
    }

    const amrex::Real src_threshold = m_src_threshold;
    for (int lev = 0; lev <= m_sim.mesh().finestLevel(); ++lev) {

        const auto& tracer_arrs = (*m_tracer)(lev).arrays();
        if (m_act_src) {
            const auto& src_arrs = (*m_act_src)(lev).const_arrays();
            amrex::ParallelFor(
                (*m_tracer)(lev), m_tracer->num_grow(),
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                    const auto src = src_arrs[nbx];
                    const amrex::Real srcmag = std::sqrt(
                        src(i, j, k, 0) * src(i, j, k, 0) +
                        src(i, j, k, 1) * src(i, j, k, 1) +
                        src(i, j, k, 2) * src(i, j, k, 2));

                    if (srcmag > src_threshold) {
                        tracer_arrs[nbx](i, j, k) = 1.0;
                    }
                });
        }

        if (m_iblank) {
            const auto& iblank_arrs = (*m_iblank)(lev).const_arrays();
            const bool tag_fringe = m_tag_fringe;
            const bool tag_hole = m_tag_hole;
            amrex::ParallelFor(
                (*m_tracer)(lev), m_tracer->num_grow(),
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                    const auto ib = iblank_arrs[nbx](i, j, k);
                    if ((tag_fringe && (ib == -1)) || (tag_hole && (ib == 0))) {
                        tracer_arrs[nbx](i, j, k) = 1.0;
                    }
                });
        }
    }
    amrex::Gpu::synchronize();
}

} // namespace amr_wind
