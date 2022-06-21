#include "amr-wind/physics/ScalarAdvection.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

ScalarAdvection::ScalarAdvection(CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{
    // Register temperature equation
    auto& teqn = sim.pde_manager().register_transport_pde("Temperature");

    // Defer getting temperature field until PDE has been registered
    m_temperature = &(teqn.fields().field);

    amrex::ParmParse pp_scalar_advection("ScalarAdvection");
    pp_scalar_advection.queryarr("location", m_loc, 0, AMREX_SPACEDIM);
    pp_scalar_advection.query("radius", m_tracer_radius);
    pp_scalar_advection.query("inner_value", m_tracer_inner);
    pp_scalar_advection.query("outer_value", m_tracer_outer);

    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.query("density", m_rho);
}

/** Initialize the velocity and temperature fields at the beginning of the
 *  simulation.
 */
void ScalarAdvection::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& density = m_density(level);
    density.setVal(m_rho);

    auto& scalars = (*m_temperature)(level);

    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();

    const amrex::Real ti = m_tracer_inner;
    const amrex::Real to = m_tracer_outer;
    const amrex::Real xc = m_loc[0];
    const amrex::Real yc = m_loc[1];
    const amrex::Real zc = m_loc[2];
    const amrex::Real radius = m_tracer_radius;

    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        auto vel = velocity.array(mfi);
        auto tracer = scalars.array(mfi);

        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
            const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
            const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
    
            vel(i, j, k, 0) = 0.0;
            vel(i, j, k, 1) = 0.0;
            vel(i, j, k, 2) = 0.0;
    
            amrex::Real r = std::sqrt(
                (x - xc) * (x - xc) + (y - yc) * (y - yc) + (z - zc) * (z - zc));
    
            if (r < radius) {
                tracer(i, j, k, 0) = ti;
            } else {
                tracer(i, j, k, 0) = to;
            }
        });
    }
}

} // namespace amr_wind
