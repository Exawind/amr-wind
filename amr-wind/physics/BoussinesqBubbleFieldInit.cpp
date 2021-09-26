#include <cmath>

#include "amr-wind/physics/BoussinesqBubbleFieldInit.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

BoussinesqBubbleFieldInit::BoussinesqBubbleFieldInit()
{
    amrex::ParmParse pp_boussinesq_bubble("BoussinesqBubble");
    pp_boussinesq_bubble.queryarr("location", m_loc, 0, AMREX_SPACEDIM);
    pp_boussinesq_bubble.query("radius", m_tracer_radius);
    pp_boussinesq_bubble.query("inner_value", m_tracer_inner);
    pp_boussinesq_bubble.query("outer_value", m_tracer_outer);

    pp_boussinesq_bubble.queryarr("init_velocity",m_init_vel,0,AMREX_SPACEDIM);

    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.query("density", m_rho);
}

void BoussinesqBubbleFieldInit::operator()(
    const amrex::Box& vbx,
    const amrex::Geometry& geom,
    const amrex::Array4<amrex::Real>& velocity,
    const amrex::Array4<amrex::Real>& density,
    const amrex::Array4<amrex::Real>& tracer,
    const amrex::Array4<amrex::Real>& mesh_fac) const
{
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();

    const amrex::Real ti = m_tracer_inner;
    const amrex::Real to = m_tracer_outer;
    const amrex::Real rho = m_rho;
    const amrex::Real xc = m_loc[0];
    const amrex::Real yc = m_loc[1];
    const amrex::Real zc = m_loc[2];
    const amrex::Real radius = m_tracer_radius;

    const amrex::Real init_v_x = m_init_vel[0];
    const amrex::Real init_v_y = m_init_vel[1];
    const amrex::Real init_v_z = m_init_vel[2];

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real x = problo[0] + (i + 0.5) * dx[0]*mesh_fac(i, j, k, 0);
        const amrex::Real y = problo[1] + (j + 0.5) * dx[1]*mesh_fac(i, j, k, 1);
        const amrex::Real z = problo[2] + (k + 0.5) * dx[2]*mesh_fac(i, j, k, 2);

        density(i, j, k) = rho;

        velocity(i, j, k, 0) = init_v_x;
        velocity(i, j, k, 1) = init_v_y;
        velocity(i, j, k, 2) = init_v_z;

        amrex::Real r = std::sqrt(
            (x - xc) * (x - xc) + (y - yc) * (y - yc) + (z - zc) * (z - zc));

        if (r < radius) {
            tracer(i, j, k, 0) = ti;
        } else {
            tracer(i, j, k, 0) = to;
        }
    });
}

} // namespace amr_wind
