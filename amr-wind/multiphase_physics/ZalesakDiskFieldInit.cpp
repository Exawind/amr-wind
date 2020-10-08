#include <cmath>

#include "amr-wind/multiphase_physics/ZalesakDiskFieldInit.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

ZalesakDiskFieldInit::ZalesakDiskFieldInit()
{

    amrex::ParmParse pp_vortex_patch("ZalesakDisk");
    pp_vortex_patch.queryarr("location", m_loc, 0, AMREX_SPACEDIM);
    pp_vortex_patch.query("radius", m_radius);
    pp_vortex_patch.query("period", m_TT);

    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.query("density", m_rho1);
}

void ZalesakDiskFieldInit::operator()(
    const amrex::Box& vbx,
    const amrex::Geometry& geom,
    const amrex::Array4<amrex::Real>& velocity,
    const amrex::Array4<amrex::Real>& density,
    const amrex::Array4<amrex::Real>& levelset) const
{

    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();

    const amrex::Real rho = m_rho1;
    const amrex::Real xc = m_loc[0];
    const amrex::Real yc = m_loc[1];
    const amrex::Real zc = m_loc[2];
    const amrex::Real radius = m_radius;

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
        const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
        const amrex::Real z = problo[1] + (k + 0.5) * dx[2];

        density(i, j, k) = rho;

        velocity(i, j, k, 0) = 2.0 * M_PI / m_TT * (0.5 - y);
        velocity(i, j, k, 1) = 2.0 * M_PI / m_TT * (x - 0.5);
        velocity(i, j, k, 2) = 0.0;

        levelset(i, j, k) =
            radius - std::sqrt(
                         (x - xc) * (x - xc) + (y - yc) * (y - yc) +
                         (z - zc) * (z - zc));
        amrex::Real d1, d2, min_dist;

        if (y < radius + yc && y > yc + radius - m_depth &&
            std::abs(x - xc) <= m_width && std::abs(z - zc) <= radius) {
            if (x > xc) {
                d1 = std::abs(xc + m_width - x);
            } else {
                d1 = std::abs(xc - m_width - x);
            }

            d2 = std::abs(y - (yc + radius - m_depth));
            min_dist = std::min(d1, d2);

            levelset(i, j, k) = -min_dist;
        }
    });
}

} // namespace amr_wind
