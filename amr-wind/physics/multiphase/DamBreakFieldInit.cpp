#include <cmath>

#include "amr-wind/physics/multiphase/DamBreakFieldInit.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

DamBreakFieldInit::DamBreakFieldInit()
{

    amrex::ParmParse pp_vortex_patch("DamBreak");
    pp_vortex_patch.queryarr("location", m_loc, 0, AMREX_SPACEDIM);
    pp_vortex_patch.query("width", m_width);
    pp_vortex_patch.query("height", m_height);
}

void DamBreakFieldInit::operator()(
    const amrex::Box& vbx,
    const amrex::Geometry& geom,
    const amrex::Array4<amrex::Real>& velocity,
    const amrex::Array4<amrex::Real>& levelset) const
{

    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();

    const amrex::Real xc = m_loc[0];
    const amrex::Real zc = m_loc[2];
    const amrex::Real width = m_width;
    const amrex::Real height = m_height;

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

        velocity(i, j, k, 0) = 0.0;
        velocity(i, j, k, 2) = 0.0;

        if (x - xc < width && z - zc < height) {
            levelset(i, j, k) = amrex::min(width - x, height - z);
        } else if (x - xc < width && z - zc > height) {
            levelset(i, j, k) = height - (z - zc);
        } else if (x - xc > width && z - zc < height) {
            levelset(i, j, k) = width - (x - xc);
        } else {
            levelset(i, j, k) = amrex::min(width - (x - xc), height - (z - zc));
        }
    });
}

} // namespace amr_wind