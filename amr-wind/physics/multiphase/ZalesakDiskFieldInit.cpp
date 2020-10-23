#include <cmath>

#include "amr-wind/physics/multiphase/ZalesakDiskFieldInit.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

ZalesakDiskFieldInit::ZalesakDiskFieldInit()
{

    amrex::ParmParse pp_vortex_patch("ZalesakDisk");
    pp_vortex_patch.queryarr("location", m_loc, 0, AMREX_SPACEDIM);
    pp_vortex_patch.query("radius", m_radius);
    pp_vortex_patch.query("period", m_TT);
}

void ZalesakDiskFieldInit::operator()(
    const amrex::Box& vbx,
    const amrex::Geometry& geom,
    const amrex::Array4<amrex::Real>& velocity,
    const amrex::Array4<amrex::Real>& levelset) const
{

    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();

    const amrex::Real xc = m_loc[0];
    const amrex::Real yc = m_loc[1];
    const amrex::Real zc = m_loc[2];
    const amrex::Real radius = m_radius;
    const amrex::Real TT = m_TT;
    const amrex::Real width = m_width;
    const amrex::Real depth = m_depth;

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
        const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

        velocity(i, j, k, 0) = 2.0 * M_PI / TT * (0.5 - y);
        velocity(i, j, k, 1) = 2.0 * M_PI / TT * (x - 0.5);
        velocity(i, j, k, 2) = 0.0;

        levelset(i, j, k) =
            radius - std::sqrt(
                         (x - xc) * (x - xc) + (y - yc) * (y - yc) +
                         (z - zc) * (z - zc));
        amrex::Real d1, d2, min_dist;

        if (y < radius + yc && y > yc + radius - depth &&
            std::abs(x - xc) <= width && std::abs(z - zc) <= radius) {
            if (x > xc) {
                d1 = std::abs(xc + width - x);
            } else {
                d1 = std::abs(xc - width - x);
            }

            d2 = std::abs(y - (yc + radius - depth));
            min_dist = amrex::min(d1, d2);

            levelset(i, j, k) = -min_dist;
        }
    });
}

} // namespace amr_wind
