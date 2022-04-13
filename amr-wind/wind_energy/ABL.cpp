#include <memory>

#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/wind_energy/ABLFieldInit.H"
#include "amr-wind/wind_energy/ABLBoundaryPlane.H"
#include "amr-wind/equation_systems/icns/source_terms/ABLForcing.H"
#include "amr-wind/equation_systems/icns/source_terms/ABLMeanBoussinesq.H"
#include "amr-wind/incflo.H"

#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"

#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldFillPatchOps.H"
#include "amr-wind/equation_systems/BCOps.H"
//#include "amr-wind/boundary_conditions/velocity_bcs.H"
#include "amr-wind/wind_energy/ShearStress.H"


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

    {
        std::string statistics_mode = "precursor";
        int dir = 2;
        amrex::ParmParse pp("ABL");
        pp.query("enable_hybrid_rl_mode", m_hybrid_rl);
        pp.query("initial_sdr_value", m_init_sdr);
        pp.query("normal_direction", dir);
        pp.query("statistics_mode", statistics_mode);
        m_stats =
            ABLStatsBase::create(statistics_mode, sim, m_abl_wall_func, dir);
    }

    // Instantiate the ABL field initializer
    m_field_init = std::make_unique<ABLFieldInit>(sim);

    // Instantiate the ABL boundary plane IO
    m_bndry_plane = std::make_unique<ABLBoundaryPlane>(sim);
}

ABL::~ABL() = default;

/** Initialize the velocity and temperature fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::ABLFieldInit
 */
void ABL::initialize_fields(int level, const amrex::Geometry& geom)
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
            level, mfi, vbx, geom, velocity.array(mfi), density.array(mfi),
            temp.array(mfi));
    }

    if (m_sim.repo().field_exists("tke")) {
        m_tke = &(m_sim.repo().get_field("tke"));
        auto& tke = (*m_tke)(level);

        for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            m_field_init->init_tke(level, mfi, vbx, geom, tke.array(mfi));
        }
    }

    if (m_sim.repo().field_exists("sdr")) {
        m_sdr = &(m_sim.repo().get_field("sdr"));
        auto& sdr = (*m_sdr)(level);

        for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            m_field_init->init_sdr(level, mfi, vbx, geom, sdr.array(mfi));
        }
    }

    if (m_sim.repo().field_exists("eps")) {
        m_eps = &(m_sim.repo().get_field("eps"));
        auto& eps = (*m_eps)(level);

        for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            m_field_init->init_eps(level, mfi, vbx, geom, eps.array(mfi));
        }
    }
}

void ABL::post_init_actions()
{
    if (m_hybrid_rl) {
        m_sdr = &(m_sim.repo().get_field("sdr"));
        m_sdr->setVal(m_init_sdr);
    }

    m_stats->post_init_actions();

    m_abl_wall_func.init_log_law_height();

    m_abl_wall_func.update_umean(
        m_stats->vel_profile(), m_stats->theta_profile());

    // Register ABL wall function for velocity
    m_velocity.register_custom_bc<ABLVelWallFunc>(m_abl_wall_func);
    (*m_temperature).register_custom_bc<ABLTempWallFunc>(m_abl_wall_func);

    m_bndry_plane->post_init_actions();
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
    const auto& vel_pa = m_stats->vel_profile();
    m_abl_wall_func.update_umean(
        m_stats->vel_profile(), m_stats->theta_profile());

    // Update the values for TKE, Omega, and Epsilon
    // if using the RANS wall model
    amrex::Orientation zlo(amrex::Direction::z, amrex::Orientation::low);
    amrex::Orientation zhi(amrex::Direction::z, amrex::Orientation::high);

    if (m_sim.repo().field_exists("tke")) {
        auto& tke = m_sim.repo().get_field("tke");
	if (tke.bc_type()[zlo] == BC::rans_wall_model ||
	    tke.bc_type()[zhi] == BC::rans_wall_model) {

	    const auto& mo = m_abl_wall_func.mo();
	    amrex::Real Cmu = m_abl_wall_func.Cmu();
	    auto tau = ShearStressAlinot(mo, Cmu);
	  
	    const amrex::Real tke_default = tau.calc_tke(); 
	    BCScalar bc_tke(tke);
	    bc_tke.update_bcvalue(tke_default);
	}
    }

    if (m_sim.repo().field_exists("sdr")) {
        auto& sdr = m_sim.repo().get_field("sdr");
	if (sdr.bc_type()[zlo] == BC::rans_wall_model ||
	    sdr.bc_type()[zhi] == BC::rans_wall_model) {

	    const auto& mo = m_abl_wall_func.mo();
	    amrex::Real Cmu = m_abl_wall_func.Cmu();
	    auto tau = ShearStressAlinot(mo, Cmu);

	    const amrex::Real sdr_default = tau.calc_omega(); 
	    BCScalar bc_sdr(sdr);
	    bc_sdr.update_bcvalue(sdr_default);
	}
    }


    if (m_abl_forcing != nullptr) {
        const amrex::Real zh = m_abl_forcing->forcing_height();
        const amrex::Real vx = vel_pa.line_average_interpolated(zh, 0);
        const amrex::Real vy = vel_pa.line_average_interpolated(zh, 1);
        // Set the mean velocities at the forcing height so that the source
        // terms can be computed during the time integration calls
        m_abl_forcing->set_mean_velocities(vx, vy);
    }

    if (m_abl_mean_bous != nullptr) {
        m_abl_mean_bous->mean_temperature_update(m_stats->theta_profile());
    }

    m_bndry_plane->pre_advance_work();
}

/** Perform tasks at the end of a new timestep
 *
 *  For ABL simulations, this method writes all plane-averaged profiles and
 *  integrated statistics to output
 */
void ABL::post_advance_work()
{
    m_stats->post_advance_work();
    m_bndry_plane->post_advance_work();
}

} // namespace amr_wind
