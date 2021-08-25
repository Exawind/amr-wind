#include "amr-wind/CFDSim.H"
#include "amr-wind/physics/TracerTag.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

TracerTag::TracerTag(CFDSim& sim)
    : m_sim(sim), m_act_src(sim.repo().get_field("actuator_src_term"))
{
    auto& teqn = sim.pde_manager().register_transport_pde("TaggingScalar");
    m_tracer = &(teqn.fields().field);

    amrex::ParmParse pp("TracerTag");
    pp.query("act_src_threshold", m_src_threshold);
}

/** Initialize the velocity and density fields at the beginning of the
 *  simulation.
 */
void TracerTag::initialize_fields(int level, const amrex::Geometry&)
{
    (*m_tracer)(level).setVal(0.0);
}

void TracerTag::post_advance_work()
{

    for (int level = 0; level <= m_sim.mesh().finestLevel(); ++level) {
        for (amrex::MFIter mfi((*m_tracer)(level)); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.tilebox();
            auto tracer = (*m_tracer)(level).array(mfi);
            auto src = m_act_src(level).array(mfi);

            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real srcmag = std::sqrt(
                        src(i, j, k, 0) * src(i, j, k, 0) +
                        src(i, j, k, 1) * src(i, j, k, 1) +
                        src(i, j, k, 2) * src(i, j, k, 2));

                    if (srcmag > m_src_threshold) tracer(i, j, k) = 1.0;
                });
        }
    }
}

} // namespace amr_wind
