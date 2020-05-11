#include "ABL.H"
#include "ABLFieldInit.H"
#include "ABLForcing.H"
#include "incflo.H"

#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "PlaneAveraging.H"

namespace amr_wind {

ABL::ABL(CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.pde_manager().icns().fields().field)
    , m_mueff(sim.pde_manager().icns().fields().mueff)
    , m_density(sim.repo().get_field("density"))
    , m_pa(2)
    , m_abl_wall_func(sim)
{
    // Register temperature equation
    // FIXME: this should be optional?
    auto& teqn = sim.pde_manager().register_transport_pde("Temperature");
    m_temperature = &(teqn.fields().field);

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
    auto& temp = (*m_temperature)(level);

    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        (*m_field_init)(
            vbx, geom, velocity.array(mfi), density.array(mfi),
            temp.array(mfi));
    }
}

void ABL::post_init_actions()
{
    m_abl_wall_func.init_log_law_height();

    // Register ABL wall function for velocity
    m_velocity.register_custom_bc<ABLVelWallFunc>(m_abl_wall_func);
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
    m_pa(
        geom, m_density.vec_ptrs(), m_velocity.vec_ptrs(), m_mueff.vec_ptrs(),
        m_temperature->vec_ptrs());

    m_abl_wall_func.update_umean(m_pa);

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
        int plot_type = 0;
        amrex::ParmParse pp("io");
        pp.query("line_plot_int", output_interval);
        pp.query("line_plot_type", plot_type);

        if ((output_interval > 0) && (time.time_index() % output_interval == 0)) {
            m_pa.plot_line(time.time_index(), time.current_time(), plot_type);
        }
    }
}

} // namespace amr_wind
