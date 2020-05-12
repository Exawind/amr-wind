#include "ABLWallFunction.H"
#include "ABL.H"
#include "tensor_ops.H"
#include "diffusion.H"

#include <cmath>

#include "AMReX_ParmParse.H"
#include "AMReX_Print.H"

namespace amr_wind {

ABLWallFunction::ABLWallFunction(const CFDSim& sim)
    : m_mesh(sim.mesh())
{
    amrex::ParmParse pp("ABL");

    pp.query("kappa", m_kappa);
    pp.query("surface_roughness_z0", m_z0);
    pp.query("normal_direction", m_direction);
    AMREX_ASSERT((0 <= m_direction) && (m_direction < AMREX_SPACEDIM));

    if (pp.contains("log_law_height")) {
        m_use_fch = false;
        pp.get("log_law_height", m_log_law_height);
    }
}

void ABLWallFunction::init_log_law_height()
{
    if (m_use_fch) {
        const auto& geom = m_mesh.Geom(0);
        m_log_law_height = (geom.ProbLo(m_direction) + 0.5 * geom.CellSize(m_direction));
    }
    {
        // fixme keeping this around to maintain perfect
        // machine zero reg tests
        // eventually turn on pre_advance_work in InitialIterations() and delete this... or make a pre_timestep_work function?
        auto geom = m_mesh.Geom();
        amrex::ParmParse pp("incflo");
        amrex::Vector<amrex::Real> vel{{0.0, 0.0, 0.0}};
        pp.queryarr("velocity", vel);
        m_umean[0] = vel[0];
        m_umean[1] = vel[1];
        m_umean[2] = vel[2];
        const amrex::Real uground = utils::vec_mag(m_umean.data());
        m_utau = m_kappa * uground / std::log(m_log_law_height / m_z0);
    }
}

void ABLWallFunction::update_umean(const PlaneAveraging& pa)
{
    m_umean[m_direction] = 0.0;
    switch (m_direction) {
    case 0:
        m_umean[1] = pa.line_velocity_ydir(m_log_law_height);
        m_umean[2] = pa.line_velocity_zdir(m_log_law_height);
        break;

    case 1:
        m_umean[0] = pa.line_velocity_xdir(m_log_law_height);
        m_umean[2] = pa.line_velocity_zdir(m_log_law_height);
        break;

    case 2:
        m_umean[0] = pa.line_velocity_xdir(m_log_law_height);
        m_umean[1] = pa.line_velocity_ydir(m_log_law_height);
        break;

    default:
        amrex::Abort("Invalid direction specified");
        break;
    }

    m_utau = m_kappa * utils::vec_mag(m_umean.data()) / (
        std::log(m_log_law_height / m_z0));
}

ABLVelWallFunc::ABLVelWallFunc(
    Field&, const ABLWallFunction& wall_func)
    : m_wall_func(wall_func)
{}

void ABLVelWallFunc::operator()(Field& velocity, const FieldState rho_state)
{
    diffusion::wall_model_bc(
        velocity, m_wall_func.utau(), utils::vec_mag(m_wall_func.umean().data()),
        rho_state);
}

} // namespace amr_wind
