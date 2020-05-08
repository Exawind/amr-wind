#include <cmath>

#include "BoussinesqBubbleFieldInit.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

BoussinesqBubbleFieldInit::BoussinesqBubbleFieldInit()
{
    amrex::ParmParse pp("BoussinesqBubble");
    pp.queryarr("location", m_loc,0,AMREX_SPACEDIM);
    pp.query("radius",m_tracer_radius);
    pp.query("inner_value",m_tracer_inner);
    pp.query("outer_value",m_tracer_outer);

    {
        amrex::ParmParse pp("incflo");
        pp.query("density", m_rho);
    }

}

void BoussinesqBubbleFieldInit::operator()(
    const amrex::Box& vbx,
    const amrex::Geometry& geom,
    const amrex::Array4<amrex::Real>& velocity,
    const amrex::Array4<amrex::Real>& density,
    const amrex::Array4<amrex::Real>& tracer) const
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

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
        const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

        density(i, j, k) = rho;

        velocity(i, j, k, 0) = 0.0;
        velocity(i, j, k, 1) = 0.0;
        velocity(i, j, k, 2) = 0.0;

        amrex::Real r = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));

        if(r < radius)
            tracer(i,j,k,0) = ti;
        else
            tracer(i,j,k,0) = to;

    });
}

} // namespace amr_wind
