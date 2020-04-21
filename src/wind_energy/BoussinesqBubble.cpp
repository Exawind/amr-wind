#include "BoussinesqBubble.H"
#include "BoussinesqBubbleFieldInit.H"
#include "incflo.H"

#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {

BoussinesqBubble::BoussinesqBubble(const SimTime& time_in, const FieldRepo& repo_in)
    : m_repo(repo_in)
{
    amrex::ParmParse pp("bb");

    pp.query("use_boussinesq", m_has_boussinesq);

    // Instantiate the BoussinesqBubble field initializer
    m_field_init.reset(new BoussinesqBubbleFieldInit());

    if (m_has_boussinesq)
        m_boussinesq.reset(new BoussinesqBuoyancy(time_in, repo_in));

}

/** Initialize the velocity and temperature fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::BoussinesqBubbleFieldInit
 */
void BoussinesqBubble::initialize_fields(
    int level,
    const amrex::Geometry& geom) const
{
    auto& velocity = m_repo.get_field("velocity")(level);
    auto& density = m_repo.get_field("density")(level);
    auto& scalars = m_repo.get_field("temperature")(level);

    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        (*m_field_init)(
            vbx, geom, velocity.array(mfi), density.array(mfi),
            scalars.array(mfi));
    }
}

void BoussinesqBubble::add_momentum_sources(
    int lev,
    FieldState fstate,
    amrex::MultiFab& vel_forces) const
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(vel_forces, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto& bx = mfi.tilebox();
        const auto& vf = vel_forces.array(mfi);

        // Boussinesq buoyancy term
        if (m_has_boussinesq) {
            (*m_boussinesq)(lev, mfi, bx, fstate, vf);
        }

    }
}

} // namespace amr_wind

