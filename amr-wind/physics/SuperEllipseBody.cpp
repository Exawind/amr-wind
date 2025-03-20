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

namespace {} // namespace

SuperEllipseBody::SuperEllipseBody(CFDSim& sim)
    : m_sim(sim)
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_body_blank(sim.repo().declare_int_field("body_blank", 1, 1, 1))
    , m_body_drag(sim.repo().declare_int_field("body_drag", 1, 1, 1))
{

        amrex::ParmParse pp(identifier());
        pp.query("body_file", m_body_file);

    m_sim.io_manager().register_output_int_var("body_drag");
    m_sim.io_manager().register_output_int_var("body_blank");

    m_body_blank.setVal(0.0);
    m_body_drag.setVal(0.0);
}

void SuperEllipseBody::initialize_fields(int level, const amrex::Geometry& geom)
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::initialize_fields");

    std::ifstream file(m_body_file, std::ios::in);
    if (file.good()) {
        // Read coordinates and orientation from body file
    }
    file.close();

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    auto& blanking = m_body_blank(level);
    auto& drag = m_body_drag(level);
    auto levelBlanking = blanking.arrays();
    auto levelDrag = drag.arrays();
    auto levelheight = body_height.arrays();
    amrex::ParallelFor(
        blanking, m_body_blank.num_grow(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
            const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
            const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

            // Determine whether point in body
            levelBlanking[nbx](i, j, k, 0) = 1;
        });
    amrex::Gpu::streamSynchronize();
    amrex::ParallelFor(
        blanking, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            if ((levelBlanking[nbx](i, j, k, 0) == 0) && (k > 0) &&
                (levelBlanking[nbx](i, j, k - 1, 0) == 1)) {
                levelDrag[nbx](i, j, k, 0) = 1;
            } else {
                levelDrag[nbx](i, j, k, 0) = 0;
            }
        });
    amrex::Gpu::streamSynchronize();
}

void SuperEllipseBody::post_init_actions()
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::post_init_actions");
    compute_body_fields();
}

void SuperEllipseBody::pre_advance_work()
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::pre_advance_work");
    compute_body_fields();
}

void SuperEllipseBody::post_regrid_actions()
{
    const int nlevels = m_sim.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        initialize_fields(lev, m_sim.repo().mesh().Geom(lev));
    }
}

void SuperEllipseBody::compute_body_fields()
{
    const int nlevels = m_sim.repo().num_active_levels();
    // Uniform, low roughness for waves
    m_bodyz0.setVal(1e-4);
    for (int level = 0; level < nlevels; ++level) {
        const auto geom = m_sim.repo().mesh().Geom(level);
        const auto& dx = geom.CellSizeArray();
        const auto& prob_lo = geom.ProbLoArray();
        auto& blanking = m_body_blank(level);
        auto& body_height = m_body_height(level);
        auto& drag = m_body_drag(level);

        auto levelBlanking = blanking.arrays();
        auto levelDrag = drag.arrays();
        auto levelHeight = body_height.arrays();

        const auto negative_wave_elevation =
            (*m_wave_negative_elevation)(level).const_arrays();
        const auto wave_vol_frac =
            (*m_wave_volume_fraction)(level).const_arrays();

        // Get body blanking from ocean waves fields
        amrex::ParallelFor(
            blanking, m_body_blank.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];                
                const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
                //Determine whether a point inside the superellipsebody
                levelBlanking[nbx](i, j, k, 0) = 1;

            });
        amrex::ParallelFor(
            blanking,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                if ((levelBlanking[nbx](i, j, k, 0) == 0) && (k > 0) &&
                    (levelBlanking[nbx](i, j, k - 1, 0) == 1)) {
                    levelDrag[nbx](i, j, k, 0) = 1;
                } else {
                    levelDrag[nbx](i, j, k, 0) = 0;
                }
            });
    }
}

} // namespace amr_wind::superellipsebody
