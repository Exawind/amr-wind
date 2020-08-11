#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/wind_energy/ABLFieldInit.H"
#include "amr-wind/equation_systems/icns/source_terms/ABLForcing.H"
#include "amr-wind/incflo.H"

#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "amr-wind/utilities/PlaneAveraging.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include "amr-wind/utilities/SecondMomentAveraging.H"
#include "amr-wind/utilities/ThirdMomentAveraging.H"

namespace amr_wind {

ABL::ABL(CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.pde_manager().icns().fields().field)
    , m_mueff(sim.pde_manager().icns().fields().mueff)
    , m_density(sim.repo().get_field("density"))
    , m_abl_wall_func(sim)
{
    // Register temperature equation
    // FIXME: this should be optional?
    auto& teqn = sim.pde_manager().register_transport_pde("Temperature");
    m_temperature = &(teqn.fields().field);

    m_stats.reset( new ABLStats(sim, m_abl_wall_func));
        
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

    if (m_field_init->add_temperature_perturbations()) {
        m_field_init->perturb_temperature(level, geom, *m_temperature);
    } else {
        temp.setVal(0.0);
    }

    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        (*m_field_init)(
            vbx, geom, velocity.array(mfi), density.array(mfi),
            temp.array(mfi));
    }

    m_stats->initialize();

    if (m_sim.repo().field_exists("tke")) {
        m_tke = &(m_sim.repo().get_field("tke"));
        auto& tke = (*m_tke)(level);
        m_field_init->init_tke(geom, tke);
    }
    

}

void ABL::post_init_actions()
{
    m_abl_wall_func.init_log_law_height();

    m_stats->calc_averages();
    
    const VelPlaneAveraging& vel_pa = m_stats->vel_plane_averaging();
    m_abl_wall_func.update_umean(vel_pa);

    // Register ABL wall function for velocity
    m_velocity.register_custom_bc<ABLVelWallFunc>(m_abl_wall_func);
}

/** Perform tasks at the beginning of a new timestep
 *
 *  For ABL simulations this method invokes the FieldPlaneAveraging class to
 *  compute spatial averages at all z-levels on the coarsest mesh (level 0).
 *
 *  The spatially averaged velocity is used to determine the current mean
 *  velocity at the forcing height (if driving pressure gradient term is active)
 *  and also determines the average friction velocity for use in the ABL wall
 *  function computation.
 */
void ABL::pre_advance_work()
{
    m_stats->calc_averages();
    const auto& vel_pa = m_stats->vel_plane_averaging();
    m_abl_wall_func.update_umean(vel_pa);

    if (m_abl_forcing != nullptr) {
        const amrex::Real zh = m_abl_forcing->forcing_height();
        const amrex::Real vx = vel_pa.line_average_interpolated(zh, 0);
        const amrex::Real vy = vel_pa.line_average_interpolated(zh, 1);
        // Set the mean velocities at the forcing height so that the source
        // terms can be computed during the time integration calls
        m_abl_forcing->set_mean_velocities(vx, vy);
    }

}

/** Perform tasks at the end of a new timestep
 *
 *  For ABL simulations, this method writes all plane-averaged profiles and
 *  integrated statistics to output 
 */
void ABL::post_advance_work()
{
    m_stats->post_advance_work();
}
    

} // namespace amr_wind
