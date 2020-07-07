#include "amr-wind/wind_energy/ABLWallFunction.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/diffusion/diffusion.H"

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
}

void ABLWallFunction::update_umean(const FieldPlaneAveraging& pa)
{
    m_umean[m_direction] = 0.0;
    switch (m_direction) {
    case 0:
        m_umean[1] = pa.line_average_interpolated(m_log_law_height, 1);
        m_umean[2] = pa.line_average_interpolated(m_log_law_height, 2);
        break;

    case 1:
        m_umean[0] = pa.line_average_interpolated(m_log_law_height, 0);
        m_umean[2] = pa.line_average_interpolated(m_log_law_height, 2);
        break;

    case 2:
        m_umean[0] = pa.line_average_interpolated(m_log_law_height, 0);
        m_umean[1] = pa.line_average_interpolated(m_log_law_height, 1);
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
