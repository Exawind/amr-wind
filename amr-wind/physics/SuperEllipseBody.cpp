#include "amr-wind/physics/SuperEllipseBody.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParReduce.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/linear_interpolation.H"

namespace amr_wind::superellipsebody {

namespace {

//! Utility function to parse inputs and return a vector instance
inline vs::Vector parse_vector(amrex::ParmParse& pp, const std::string& key)
{
    amrex::Vector<amrex::Real> tmp;
    pp.getarr(key.data(), tmp);
    AMREX_ALWAYS_ASSERT(tmp.size() == 3U);

    return vs::Vector{tmp[0], tmp[1], tmp[2]};
}


} // namespace

SuperEllipseBody::SuperEllipseBody(CFDSim& sim)
    : m_sim(sim)
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_body_blank(sim.repo().declare_int_field("terrain_blank", 1, 1, 1))
{

        amrex::ParmParse pp(identifier());
        pp.query("body_file", m_body_file);
        m_dim = parse_vector(pp, "dimensions");

    m_sim.io_manager().register_output_int_var("body_drag");

    m_body_blank.setVal(0.0);
}

void SuperEllipseBody::initialize_fields(int level, const amrex::Geometry& geom)
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::initialize_fields");

    std::ifstream file(m_body_file, std::ios::in);
    if (file.good()) {
        // Read coordinates and orientation from body file
        // Just read 12 numbers one after the other. First 3 represent location, next 9 represent orientation
        file >> m_loc[0] >> m_loc[1] >> m_loc[2];
        file >> m_orient.xx() >> m_orient.xy() >> m_orient.xz() >> m_orient.yx()
             >> m_orient.yy() >> m_orient.yz() >> m_orient.zx() >> m_orient.zy()
             >> m_orient.zz();
    }
    file.close();

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    auto& blanking = m_body_blank(level);
    auto levelBlanking = blanking.arrays();

    vs::Vector dim = m_dim;
    vs::Tensor orient = m_orient;
    vs::Vector loc = m_loc;

    amrex::ParallelFor(
        blanking, m_body_blank.num_grow(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
            const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
            const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

            // Transform to orientation inside 
            vs::Vector xtransform = orient & (vs::Vector(x, y, z) - loc);

            // Determine whether point in body
            levelBlanking[nbx](i, j, k, 0) = static_cast<int>( (xtransform.x()/dim.x())^6 + (xtransform.y()/dim.y())^6 + (xtransform.z()/dim.z())^6 - 1.0 );
        });
    amrex::Gpu::streamSynchronize();
}

void SuperEllipseBody::post_regrid_actions()
{

    const int nlevels = m_sim.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        initialize_fields(lev, m_sim.repo().mesh().Geom(lev));
    }
}

} // namespace amr_wind::superellipsebody
