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
    , m_bodyz0(sim.repo().declare_field("bodyz0", 1, 1, 1))
    , m_body_height(sim.repo().declare_field("body_height", 1, 1, 1))
{

    m_body_is_waves = sim.physics_manager().contains("OceanWaves") &&
                         !sim.repo().field_exists("vof");

    if (!m_body_is_waves) {
        amrex::ParmParse pp(identifier());
        pp.query("body_file", m_body_file);
        pp.query("roughness_file", m_roughness_file);
    } else {
        m_wave_volume_fraction = &m_repo.get_field(m_wave_volume_fraction_name);
        m_wave_negative_elevation =
            &m_repo.get_field(m_wave_negative_elevation_name);
    }

    m_sim.io_manager().register_output_int_var("body_drag");
    m_sim.io_manager().register_output_int_var("body_blank");
    m_sim.io_manager().register_io_var("bodyz0");
    m_sim.io_manager().register_io_var("body_height");

    m_body_blank.setVal(0.0);
    m_body_drag.setVal(0.0);
    m_bodyz0.set_default_fillpatch_bc(m_sim.time());
    m_body_height.set_default_fillpatch_bc(m_sim.time());
}

void SuperEllipseBody::initialize_fields(int level, const amrex::Geometry& geom)
{
    if (m_body_is_waves) {
        return;
    }

    BL_PROFILE("amr-wind::" + this->identifier() + "::initialize_fields");

    //! Reading the Body Coordinates from  file
    amrex::Vector<amrex::Real> xbody;
    amrex::Vector<amrex::Real> ybody;
    amrex::Vector<amrex::Real> zbody;
    ioutils::read_flat_grid_file(m_body_file, xbody, ybody, zbody);

    // No checks for the file as it is optional currently
    amrex::Vector<amrex::Real> xrough;
    amrex::Vector<amrex::Real> yrough;
    amrex::Vector<amrex::Real> z0rough;
    std::ifstream file(m_roughness_file, std::ios::in);
    if (file.good()) {
        ioutils::read_flat_grid_file(m_roughness_file, xrough, yrough, z0rough);
    }
    file.close();

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    auto& blanking = m_body_blank(level);
    auto& bodyz0 = m_bodyz0(level);
    auto& body_height = m_body_height(level);
    auto& drag = m_body_drag(level);
    const auto xbody_size = xbody.size();
    const auto ybody_size = ybody.size();
    const auto zbody_size = zbody.size();
    amrex::Gpu::DeviceVector<amrex::Real> d_xbody(xbody_size);
    amrex::Gpu::DeviceVector<amrex::Real> d_ybody(ybody_size);
    amrex::Gpu::DeviceVector<amrex::Real> d_zbody(zbody_size);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, xbody.begin(), xbody.end(),
        d_xbody.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, ybody.begin(), ybody.end(),
        d_ybody.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, zbody.begin(), zbody.end(),
        d_zbody.begin());
    const auto* xbody_ptr = d_xbody.data();
    const auto* ybody_ptr = d_ybody.data();
    const auto* zbody_ptr = d_zbody.data();
    // Copy Roughness to gpu
    const auto xrough_size = xrough.size();
    const auto yrough_size = yrough.size();
    const auto z0rough_size = z0rough.size();
    amrex::Gpu::DeviceVector<amrex::Real> d_xrough(xrough_size);
    amrex::Gpu::DeviceVector<amrex::Real> d_yrough(yrough_size);
    amrex::Gpu::DeviceVector<amrex::Real> d_z0rough(z0rough_size);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, xrough.begin(), xrough.end(),
        d_xrough.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, yrough.begin(), yrough.end(),
        d_yrough.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, z0rough.begin(), z0rough.end(),
        d_z0rough.begin());
    const auto* xrough_ptr = d_xrough.data();
    const auto* yrough_ptr = d_yrough.data();
    const auto* z0rough_ptr = d_z0rough.data();
    auto levelBlanking = blanking.arrays();
    auto levelDrag = drag.arrays();
    auto levelz0 = bodyz0.arrays();
    auto levelheight = body_height.arrays();
    amrex::ParallelFor(
        blanking, m_body_blank.num_grow(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
            const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
            const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
            const amrex::Real bodyHt = interp::bilinear(
                xbody_ptr, xbody_ptr + xbody_size, ybody_ptr,
                ybody_ptr + ybody_size, zbody_ptr, x, y);
            levelBlanking[nbx](i, j, k, 0) =
                static_cast<int>((z <= bodyHt) && (z > prob_lo[2]));
            levelheight[nbx](i, j, k, 0) = bodyHt;

            amrex::Real roughz0 = 0.1;
            if (xrough_size > 0) {
                roughz0 = interp::bilinear(
                    xrough_ptr, xrough_ptr + xrough_size, yrough_ptr,
                    yrough_ptr + yrough_size, z0rough_ptr, x, y);
            }
            levelz0[nbx](i, j, k, 0) = roughz0;
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
    if (!m_body_is_waves) {
        return;
    }
    BL_PROFILE("amr-wind::" + this->identifier() + "::post_init_actions");
    convert_waves_to_body_fields();
}

void SuperEllipseBody::pre_advance_work()
{
    if (!m_body_is_waves) {
        return;
    }
    BL_PROFILE("amr-wind::" + this->identifier() + "::pre_advance_work");
    convert_waves_to_body_fields();
}

void SuperEllipseBody::post_regrid_actions()
{
    if (m_body_is_waves) {
        convert_waves_to_body_fields();
    } else {
        const int nlevels = m_sim.repo().num_active_levels();
        for (int lev = 0; lev < nlevels; ++lev) {
            initialize_fields(lev, m_sim.repo().mesh().Geom(lev));
        }
    }
}

void SuperEllipseBody::convert_waves_to_body_fields()
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
                const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
                levelBlanking[nbx](i, j, k, 0) = static_cast<int>(
                    (wave_vol_frac[nbx](i, j, k) >= 0.5) && (z > prob_lo[2]));
                levelHeight[nbx](i, j, k, 0) =
                    -negative_wave_elevation[nbx](i, j, k);
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
