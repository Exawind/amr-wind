#include "ABL.H"
#include "ABLFieldInit.H"
#include "incflo.H"

#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "PlaneAveraging.H"

namespace amr_wind {

ABL::ABL(const SimTime& time, incflo* incflo_in)
    : m_time(time), m_incflo(incflo_in)
{
    amrex::ParmParse pp("abl");

    pp.query("use_boussinesq", m_has_boussinesq);
    pp.query("coriolis_effect", m_has_coriolis);
    pp.query("abl_forcing", m_has_driving_dpdx);

    // Instantiate the ABL field initializer
    m_field_init.reset(new ABLFieldInit());
    m_abl_wall_func.reset(new ABLWallFunction());

    if (m_has_boussinesq)
        m_boussinesq.reset(new BoussinesqBuoyancy());

    if (m_has_driving_dpdx)
        m_abl_forcing.reset(new ABLForcing(m_time));

    if (m_has_coriolis)
        m_coriolis.reset(new CoriolisForcing());
}

/** Initialize the velocity and temperature fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::ABLFieldInit
 */
void ABL::initialize_fields(
    const amrex::Geometry& geom,
    LevelData& leveldata) const
{
    auto& velocity = leveldata.velocity;
    auto& density = leveldata.density;
    auto& scalars = leveldata.tracer;

    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        (*m_field_init)(
            vbx, geom, velocity.array(mfi), density.array(mfi),
            scalars.array(mfi));
    }
}

/** Add ABL-related source terms to the momentum equation
 */
void ABL::add_momentum_sources(
    const amrex::Geometry& /* geom */,
    const LevelData& leveldata,
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
            const auto& scalars_old = leveldata.tracer_o.const_array(mfi);
            const auto& scalars_new = leveldata.tracer.const_array(mfi);
            (*m_boussinesq)(bx, scalars_old, scalars_new, vf);
        }

        // Coriolis term
        if (m_has_coriolis) {
            const auto& vel = leveldata.velocity.const_array(mfi);
            (*m_coriolis)(bx, vel, vf);
        }

        // Driving pressure gradient term
        if (m_has_driving_dpdx) (*m_abl_forcing)(bx, vf);
    }
}

void ABL::add_momentum_sources(
    const amrex::Geometry& /* geom */,
    const amrex::MultiFab& /* density */,
    const amrex::MultiFab& velocity,
    const amrex::MultiFab& scalars_old,
    const amrex::MultiFab& scalars_new,
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
            const auto& scal_o = scalars_old.const_array(mfi);
            const auto& scal_n = scalars_new.const_array(mfi);
            (*m_boussinesq)(bx, scal_o, scal_n, vf);
        }

        // Coriolis term
        if (m_has_coriolis) {
            const auto& vel = velocity.const_array(mfi);
            (*m_coriolis)(bx, vel, vf);
        }

        // Driving pressure gradient term
        if (m_has_driving_dpdx) (*m_abl_forcing)(bx, vf);
    }
}


/** Perform tasks at the beginning of a new timestep
 *
 *  For ABL simulations this method invokes the PlaneAveraging class to
 *  compute spatial averages at all z-levels on the coarsest mesh (level 0).
 *
 *  The spatially averaged velocity is used to determine the current mean
 *  velocity at the forcing height (if driving pressure gradient term is active)
 *  and also determines the average friction velocity for use in the ABL wall
 *  function computation.
 */
void ABL::pre_advance_work()
{
    // Spatial averaging on z-planes
    constexpr int direction = 2;
    const auto& geom = m_incflo->Geom(0);

    // TODO: Promote this to a class member
    PlaneAveraging pa(geom, *m_incflo->leveldata_vec()[0], direction);
    {
        // First cell height
        const amrex::Real fch = geom.ProbLo(direction) + 0.5 * geom.CellSize(direction);
        const amrex::Real vx = pa.line_velocity_xdir(fch);
        const amrex::Real vy = pa.line_velocity_ydir(fch);

        const amrex::Real uground = std::sqrt(vx * vx + vy * vy);
        const amrex::Real utau = m_abl_wall_func->utau(uground, fch);
        // TODO: For now retain wall function logic in incflo
        m_incflo->set_abl_friction_vels(uground, utau);
    }

    if (m_has_driving_dpdx) {
        const amrex::Real zh = m_abl_forcing->forcing_height();
        const amrex::Real vx = pa.line_velocity_xdir(zh);
        const amrex::Real vy = pa.line_velocity_ydir(zh);
        // Set the mean velocities at the forcing height so that the source
        // terms can be computed during the time integration calls
        m_abl_forcing->set_mean_velocities(vx, vy);
    }

    {
        // TODO: This should be handled by PlaneAveraging
        int output_interval = 1;
        amrex::ParmParse pp("amr");
        pp.query("line_plot_int", output_interval);

        if ((output_interval > 0) && (m_time.time_index() % output_interval == 0)) {
            pa.plot_line_text("line_plot.txt", m_time.time_index(), m_time.current_time());
        }
    }
}

} // namespace amr_wind
