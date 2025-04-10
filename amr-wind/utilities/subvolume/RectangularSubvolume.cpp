#include "amr-wind/utilities/subvolume/RectangularSubvolume.H"
#include "amr-wind/utilities/constants.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::subvolume {

RectangularSubvolume::RectangularSubvolume(const CFDSim& sim) : m_sim(sim) {}

RectangularSubvolume::~RectangularSubvolume() = default;

void RectangularSubvolume::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);

    pp.getarr("origin", m_origin);
    pp.getarr("num_points", m_npts_vec);
    pp.getarr("dx_vec", m_dx_vec);

    m_chunk_size_vec.resize(AMREX_SPACEDIM);
    m_chunk_size_vec[0] = m_sim.mesh().maxGridSize(0)[0];
    m_chunk_size_vec[1] = m_sim.mesh().maxGridSize(0)[1];
    m_chunk_size_vec[2] = m_sim.mesh().maxGridSize(0)[2];
    pp.queryarr("chunk_size", m_chunk_size_vec);
}

void RectangularSubvolume::evaluate_inputs()
{
    bool found = false;
    const auto& geom = m_sim.mesh().Geom();
    for (int i = 0; i < m_sim.repo().num_active_levels(); i++) {
        if (!found) {
            if (std::abs(m_dx_vec[0] - geom[i].CellSize(0)) < 1e-4 &&
                std::abs(m_dx_vec[1] - geom[i].CellSize(1)) < 1e-4 &&
                std::abs(m_dx_vec[2] - geom[i].CellSize(2)) < 1e-4) {

                amrex::Print()
                    << "RectangularSubvolume " + m_label +
                           ": Resolution specified matches that of level "
                    << i << std::endl;
                found = true;
                m_lev_for_sub = i;
            }
        }
    }

    if (!found) {
        amrex::Abort(
            "RectangularSubvolume " + m_label +
            ": Resolution specified for subvolume does not match the "
            "resolution "
            "of any of the mesh levels.");
    }

    // **************************************************************
    // Now that we know which level we're at, we can figure out which (i,j,k)
    // the origin corresponds to Note we use 1.0001 as a fudge factor since the
    // division of two reals --> integer will do a floor
    // **************************************************************
    int i0 = static_cast<int>(
        (m_origin[0] - geom[m_lev_for_sub].ProbLo(0)) * 1.0001 /
        geom[m_lev_for_sub].CellSize(0));
    int j0 = static_cast<int>(
        (m_origin[1] - geom[m_lev_for_sub].ProbLo(1)) * 1.0001 /
        geom[m_lev_for_sub].CellSize(1));
    int k0 = static_cast<int>(
        (m_origin[2] - geom[m_lev_for_sub].ProbLo(2)) * 1.0001 /
        geom[m_lev_for_sub].CellSize(2));

    found = false;
    if (std::abs(
            geom[m_lev_for_sub].ProbLo(0) +
            i0 * geom[m_lev_for_sub].CellSize(0) - m_origin[0]) < 1e-4 &&
        std::abs(
            geom[m_lev_for_sub].ProbLo(1) +
            j0 * geom[m_lev_for_sub].CellSize(1) - m_origin[1]) < 1e-4 &&
        std::abs(
            geom[m_lev_for_sub].ProbLo(2) +
            k0 * geom[m_lev_for_sub].CellSize(2) - m_origin[2]) < 1e-4) {
        amrex::Print()
            << "RectangularSubvolume " + m_label +
                   ": Specified origin is the lower left corner of cell "
            << amrex::IntVect(i0, j0, k0) << std::endl;
        found = true;
    }

    if (!found) {
        amrex::Abort(
            "RectangularSubvolume " + m_label +
            ": Origin specified does not correspond to a node at this level.");
    }

    amrex::Box domain(geom[m_lev_for_sub].Domain());

    amrex::Box bx(
        amrex::IntVect(i0, j0, k0),
        amrex::IntVect(
            i0 + m_npts_vec[0] - 1, j0 + m_npts_vec[1] - 1,
            k0 + m_npts_vec[2] - 1));
    amrex::Print() << "RectangularSubvolume " + m_label + ": Box requested is "
                   << bx << std::endl;

    if (!domain.contains(bx)) {
        amrex::Abort(
            "RectangularSubvolume " + m_label +
            ": Box requested is larger than the existing domain");
    }

    amrex::IntVect chunk_size(
        m_chunk_size_vec[0], m_chunk_size_vec[1], m_chunk_size_vec[2]);

    amrex::BoxArray ba(bx);
    ba.maxSize(chunk_size);

    amrex::Print() << "RectangularSubvolume " + m_label + ": BoxArray is " << ba
                   << std::endl;

    m_ba = ba;
}

} // namespace amr_wind::subvolume
