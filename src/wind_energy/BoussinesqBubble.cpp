#include "BoussinesqBubble.H"
#include "BoussinesqBubbleFieldInit.H"
#include "incflo.H"

#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {

BoussinesqBubbleOld::BoussinesqBubbleOld(incflo* incflo_in)
    : m_incflo(incflo_in)
{
    amrex::ParmParse pp("bb");

    pp.query("use_boussinesq", m_has_boussinesq);

    // Instantiate the BoussinesqBubble field initializer
    m_field_init.reset(new BoussinesqBubbleFieldInit());

    if (m_has_boussinesq)
        m_boussinesq.reset(new BoussinesqBuoyancyOld());

}

/** Initialize the velocity and temperature fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::BoussinesqBubbleFieldInit
 */
void BoussinesqBubbleOld::initialize_fields(
    int level,
    const amrex::Geometry& geom) const
{
    auto& velocity = m_incflo->velocity()(level);
    auto& density = m_incflo->density()(level);
    auto& scalars = m_incflo->tracer()(level);

    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        (*m_field_init)(
            vbx, geom, velocity.array(mfi), density.array(mfi),
            scalars.array(mfi));
    }
}

void BoussinesqBubbleOld::add_momentum_sources(
    const amrex::Geometry& /* geom */,
    const amrex::MultiFab& /* density */,
    const amrex::MultiFab& /* velocity */,
    const amrex::MultiFab& scalars,
    amrex::MultiFab& vel_forces) const
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(vel_forces, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const auto& bx = mfi.tilebox();
        const auto& vf = vel_forces.array(mfi);

        // Boussinesq buoyancy term
        if (m_has_boussinesq) {
            const auto& scal = scalars.const_array(mfi);
            (*m_boussinesq)(bx, scal, vf);
        }
    }
}

BoussinesqBubble::BoussinesqBubble(const CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
    , m_temperature(sim.repo().get_field("temperature"))
{
    // Instantiate the BoussinesqBubble field initializer
    m_field_init.reset(new BoussinesqBubbleFieldInit());
}

/** Initialize the velocity and temperature fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::BoussinesqBubbleFieldInit
 */
void BoussinesqBubble::initialize_fields(
    int level,
    const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& density = m_density(level);
    auto& scalars = m_temperature(level);

    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        (*m_field_init)(
            vbx, geom, velocity.array(mfi), density.array(mfi),
            scalars.array(mfi));
    }
}

} // namespace amr_wind
