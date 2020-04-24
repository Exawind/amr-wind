#include "ABL.H"
#include "ABLFieldInit.H"
#include "incflo.H"

#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "PlaneAveraging.H"

namespace amr_wind {

ABLOld::ABLOld(const SimTime& time, incflo* incflo_in)
    : m_time(time), m_incflo(incflo_in), m_pa(2)
{
    amrex::ParmParse pp("abl");

    pp.query("use_boussinesq", m_has_boussinesq);
    pp.query("coriolis_effect", m_has_coriolis);
    pp.query("abl_forcing", m_has_driving_dpdx);

    // Instantiate the ABL field initializer
    m_field_init.reset(new ABLFieldInit());
    m_abl_wall_func.reset(new ABLWallFunction());

    if (m_has_boussinesq)
        m_boussinesq.reset(new BoussinesqBuoyancyOld());

    if (m_has_driving_dpdx)
        m_abl_forcing.reset(new ABLForcingOld(m_time));

    if (m_has_coriolis)
        m_coriolis.reset(new CoriolisForcingOld());

    {
        // fixme keeping this around to maintain perfect
        // machine zero reg tests
        // eventually turn on pre_advance_work in InitialIterations() and delete this... or make a pre_timestep_work function?
        constexpr int direction = 2;
        auto geom = m_incflo->Geom();
        // First cell height
        const amrex::Real fch = geom[0].ProbLo(direction) + 0.5 * geom[0].CellSize(direction);
        amrex::Real vx = 0.0;
        amrex::Real vy = 0.0;
        amrex::ParmParse pp("incflo");
        pp.query("ic_u", vx);
        pp.query("ic_v", vy);
        const amrex::Real uground = std::sqrt(vx * vx + vy * vy);
        const amrex::Real utau = m_abl_wall_func->utau(uground, fch);
        m_incflo->set_abl_friction_vels(uground, utau);
    }
}

/** Initialize the velocity and temperature fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::ABLFieldInit
 */
void ABLOld::initialize_fields(
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

void ABLOld::add_momentum_sources(
    const amrex::Geometry& /* geom */,
    const amrex::MultiFab& /* density */,
    const amrex::MultiFab& ,
    const amrex::MultiFab& ,
    amrex::MultiFab& ) const
{
#if 0
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

        // Coriolis term
        if (m_has_coriolis) {
            const auto& vel = velocity.const_array(mfi);
            (*m_coriolis)(bx, vel, vf);
        }

        // Driving pressure gradient term
        if (m_has_driving_dpdx) (*m_abl_forcing)(bx, vf);
    }
#endif
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
void ABLOld::pre_advance_work()
{
    BL_PROFILE("amr-wind::ABL::pre_advance_work")
    // Spatial averaging on z-planes
    constexpr int direction = 2;
    auto geom = m_incflo->Geom();

    m_pa(geom, m_incflo->velocity().vec_ptrs(), m_incflo->tracer().vec_ptrs());

    {
        // First cell height
        const amrex::Real fch = geom[0].ProbLo(direction) + 0.5 * geom[0].CellSize(direction);
        const amrex::Real vx = m_pa.line_velocity_xdir(fch);
        const amrex::Real vy = m_pa.line_velocity_ydir(fch);

        const amrex::Real uground = std::sqrt(vx * vx + vy * vy);
        const amrex::Real utau = m_abl_wall_func->utau(uground, fch);
        // TODO: For now retain wall function logic in incflo
        m_incflo->set_abl_friction_vels(uground, utau);
    }

    if (m_has_driving_dpdx) {
        const amrex::Real zh = m_abl_forcing->forcing_height();
        const amrex::Real vx = m_pa.line_velocity_xdir(zh);
        const amrex::Real vy = m_pa.line_velocity_ydir(zh);
        // Set the mean velocities at the forcing height so that the source
        // terms can be computed during the time integration calls
        m_abl_forcing->set_mean_velocities(vx, vy);
    }

#if 0
    {
        // TODO: This should be handled by PlaneAveraging
        int output_interval = 1;
        amrex::ParmParse pp("amr");
        pp.query("line_plot_int", output_interval);

        if ((output_interval > 0) && (m_time.time_index() % output_interval == 0)) {
            m_pa.plot_line_text("line_plot.txt", m_time.time_index(), m_time.current_time());
        }
    }
#endif
}

ABL::ABL(const CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
    , m_temperature(sim.repo().get_field("temperature"))
    , m_pa(2)
{
    // Instantiate the ABL field initializer
    m_field_init.reset(new ABLFieldInit());
}

ABL::~ABL() = default;

/** Initialize the velocity and temperature fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::ABLFieldInit
 */
void ABL::initialize_fields(
    int level,
    const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& density = m_density(level);
    auto& temp = m_temperature(level);

    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        (*m_field_init)(
            vbx, geom, velocity.array(mfi), density.array(mfi),
            temp.array(mfi));
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
    const auto& time = m_sim.time();
    const auto& geom = m_sim.mesh().Geom();
    m_pa(geom, m_velocity.vec_ptrs(), m_temperature.vec_ptrs());

    if (m_abl_forcing != nullptr) {
        const amrex::Real zh = m_abl_forcing->forcing_height();
        const amrex::Real vx = m_pa.line_velocity_xdir(zh);
        const amrex::Real vy = m_pa.line_velocity_ydir(zh);
        // Set the mean velocities at the forcing height so that the source
        // terms can be computed during the time integration calls
        m_abl_forcing->set_mean_velocities(vx, vy);
    }

    {
        // TODO: This should be handled by PlaneAveraging
        int output_interval = 1;
        amrex::ParmParse pp("amr");
        pp.query("line_plot_int", output_interval);

        if ((output_interval > 0) && (time.time_index() % output_interval == 0)) {
            m_pa.plot_line_text("line_plot.txt", time.time_index(), time.current_time());
        }
    }
}

} // namespace amr_wind
