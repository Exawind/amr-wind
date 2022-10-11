#include "amr-wind/utilities/tagging/CylinderRefiner.H"
#include "AMReX_ParmParse.H"

// Adapted from OpenFOAM/src/meshTools/sets/cellSources/cylinderAnnulusToCell

namespace amr_wind::tagging {

CylinderRefiner::CylinderRefiner(
    const CFDSim& /*unused*/, const std::string& key)
{
    amrex::ParmParse pp(key);

    amrex::Vector<amrex::Real> tmp_vec;

    // Start of the cylinder
    pp.getarr("start", tmp_vec);
    AMREX_ALWAYS_ASSERT(tmp_vec.size() == AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        m_start[i] = tmp_vec[i];
    }

    // End point of cylinder
    pp.getarr("end", tmp_vec);
    AMREX_ALWAYS_ASSERT(tmp_vec.size() == AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        m_end[i] = tmp_vec[i];
    }

    // Outer radial extent of the cylinder, always read in from input file
    pp.get("outer_radius", m_outer_radius);
    // Optional inner radius to restrict tagging to an annulus of the cylinder
    pp.query("inner_radius", m_inner_radius);
}

void CylinderRefiner::operator()(
    const amrex::Box& bx,
    const amrex::Geometry& geom,
    const amrex::Array4<amrex::TagBox::TagType>& tag) const
{
    const auto axis = m_end - m_start;
    const auto start = m_start;
    const auto outer = m_outer_radius * m_outer_radius;
    const auto inner = m_inner_radius * m_inner_radius;
    const auto magax = vs::mag_sqr(axis);
    const auto& problo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
        const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

        // Position vector of the cell center
        const vs::Vector pt(x, y, z);
        // Vector relative to the origin of the cylinder
        const vs::Vector pvec = pt - start;
        // Component along the cylinder axis
        const amrex::Real daxis = pvec & axis;

        // Check if the point lies in between the cylinder extents along axis
        if ((daxis >= 0) && (daxis <= magax)) {
            const amrex::Real d2 = (pvec & pvec) - (daxis * daxis) / magax;

            // Check if the cell center lies within the radius specified
            if ((d2 <= outer) && (d2 >= inner)) {
                tag(i, j, k) = amrex::TagBox::SET;
            }
        }
    });
}

} // namespace amr_wind::tagging
