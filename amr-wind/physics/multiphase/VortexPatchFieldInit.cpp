#include <cmath>

#include "amr-wind/physics/multiphase/VortexPatchFieldInit.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

VortexPatchFieldInit::VortexPatchFieldInit()
{

    amrex::ParmParse pp_vortex_patch("VortexPatch");
    pp_vortex_patch.queryarr("location", m_loc, 0, AMREX_SPACEDIM);
    pp_vortex_patch.query("radius", m_radius);
}

void VortexPatchFieldInit::operator()(
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

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
        const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

        velocity(i, j, k, 0) = 2.0 * std::sin(M_PI * x) * std::sin(M_PI * x) *
                               std::sin(2.0 * M_PI * y) *
                               std::sin(2.0 * M_PI * z);
        velocity(i, j, k, 1) = -std::sin(M_PI * y) * std::sin(M_PI * y) *
                               std::sin(2.0 * M_PI * x) *
                               std::sin(2.0 * M_PI * z);
        velocity(i, j, k, 2) = -std::sin(M_PI * z) * std::sin(M_PI * z) *
                               std::sin(2.0 * M_PI * x) *
                               std::sin(2.0 * M_PI * y);

        levelset(i, j, k) =
            radius - std::sqrt(
                         (x - xc) * (x - xc) + (y - yc) * (y - yc) +
                         (z - zc) * (z - zc));
    });
}

} // namespace amr_wind
