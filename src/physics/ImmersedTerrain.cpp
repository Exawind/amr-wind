#include "src/physics/ImmersedTerrain.H"
#include "src/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "src/utilities/IOManager.H"
#include "src/utilities/io_utils.H"
#include "src/utilities/linear_interpolation.H"
#include "AMReX_REAL.H"

#include <fstream>

using namespace amrex::literals;

namespace kynema_sgf::immersedterrain {

namespace {
//! Fractions within this tolerance of 0 or 1 are snapped, so that the smooth
//! distance-function blanking still yields exactly fluid / solid cells away
//! from the surface and the mask search has a clean threshold to work with.
constexpr amrex::Real fraction_tol = 1.0e-3_rt;
} // namespace

ImmersedTerrain::ImmersedTerrain(CFDSim& sim)
    : m_sim(sim)
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_terrain_fraction(sim.repo().declare_field("terrain_fraction", 1, 1, 1))
    , m_terrain_mask(sim.repo().declare_int_field("terrain_mask", 1, 1, 1))
    , m_terrain_surface(
          sim.repo().declare_field("terrain_surface", AMREX_SPACEDIM, 1, 1))
    , m_terrain_roughness(
          sim.repo().declare_field("terrain_roughness", 1, 1, 1))
{
    amrex::ParmParse pp(identifier());
    pp.query("terrain_file", m_terrain_file);
    pp.query("roughness_file", m_roughness_file);
    pp.query("uniform_roughness", m_uniform_z0);
    if (pp.contains("uniform_roughness") && pp.contains("roughness_file")) {
        amrex::Print()
            << "Warning: Both uniform_roughness and roughness_file are "
               "specified. Roughness file values will override uniform "
               "roughness provided.\n";
    }
    pp.query("blanking_method", m_blanking_method);
    pp.query("smoothing_length", m_smoothing_length);
    pp.query("solid_threshold", m_solid_threshold);
    if (m_blanking_method != "volume_fraction" &&
        m_blanking_method != "distance_function") {
        amrex::Abort(
            identifier() +
            ".blanking_method must be volume_fraction or "
            "distance_function, got " +
            m_blanking_method);
    }

    pp.query("implicit_projection", m_implicit_projection);
    if (m_implicit_projection) {
        // Same coefficient the explicit source term would use
        amrex::ParmParse pp_drag("ImmersedDragForcing");
        pp_drag.query("drag_coefficient", m_drag_coefficient);
        m_terrain_drag_rate =
            &sim.repo().declare_field("terrain_drag_rate", 1, 1, 1);
        m_terrain_drag_rate->setVal(0.0_rt);
        m_terrain_drag_rate->set_default_fillpatch_bc(m_sim.time());
        m_sim.io_manager().register_io_var("terrain_drag_rate");
        amrex::Print() << identifier()
                       << ": immersed drag applied implicitly in the "
                          "projections with C_d = "
                       << m_drag_coefficient << "\n";
    }

    m_sim.io_manager().register_output_int_var("terrain_mask");
    m_sim.io_manager().register_io_var("terrain_fraction");
    m_sim.io_manager().register_io_var("terrain_surface");
    m_sim.io_manager().register_io_var("terrain_roughness");

    m_terrain_fraction.setVal(0.0_rt);
    m_terrain_mask.setVal(mask_fluid);
    m_terrain_surface.setVal(0.0_rt);
    m_terrain_roughness.setVal(m_uniform_z0);
    m_terrain_fraction.set_default_fillpatch_bc(m_sim.time());
    m_terrain_surface.set_default_fillpatch_bc(m_sim.time());
    m_terrain_roughness.set_default_fillpatch_bc(m_sim.time());
}

void ImmersedTerrain::initialize_fields(int level, const amrex::Geometry& geom)
{
    BL_PROFILE("kynema-sgf::" + this->identifier() + "::initialize_fields");

    //! Terrain coordinates from file
    amrex::Vector<amrex::Real> xterrain;
    amrex::Vector<amrex::Real> yterrain;
    amrex::Vector<amrex::Real> zterrain;
    ioutils::read_flat_grid_file(m_terrain_file, xterrain, yterrain, zterrain);

    //! Roughness file is optional
    amrex::Vector<amrex::Real> xrough;
    amrex::Vector<amrex::Real> yrough;
    amrex::Vector<amrex::Real> z0rough;
    {
        std::ifstream file(m_roughness_file, std::ios::in);
        if (file.good()) {
            ioutils::read_flat_grid_file(
                m_roughness_file, xrough, yrough, z0rough);
        }
    }

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    auto& fraction = m_terrain_fraction(level);
    auto& mask = m_terrain_mask(level);
    auto& surface = m_terrain_surface(level);
    auto& roughness = m_terrain_roughness(level);

    // Copy terrain to device
    const auto xterrain_size = xterrain.size();
    const auto yterrain_size = yterrain.size();
    const auto zterrain_size = zterrain.size();
    amrex::Gpu::DeviceVector<amrex::Real> d_xterrain(xterrain_size);
    amrex::Gpu::DeviceVector<amrex::Real> d_yterrain(yterrain_size);
    amrex::Gpu::DeviceVector<amrex::Real> d_zterrain(zterrain_size);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, xterrain.begin(), xterrain.end(),
        d_xterrain.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, yterrain.begin(), yterrain.end(),
        d_yterrain.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, zterrain.begin(), zterrain.end(),
        d_zterrain.begin());
    const auto* xterrain_ptr = d_xterrain.data();
    const auto* yterrain_ptr = d_yterrain.data();
    const auto* zterrain_ptr = d_zterrain.data();

    // Copy roughness to device
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

    auto frac_arrs = fraction.arrays();
    auto mask_arrs = mask.arrays();
    const bool has_rate = (m_terrain_drag_rate != nullptr);
    auto rate_arrs = has_rate ? (*m_terrain_drag_rate)(level).arrays()
                              : amrex::MultiArray4<amrex::Real>();
    const amrex::Real drag_rate_solid = m_drag_coefficient / dx[2];
    auto surf_arrs = surface.arrays();
    auto z0_arrs = roughness.arrays();

    const amrex::Real uniform_z0 = m_uniform_z0;
    const bool use_distance_function =
        (m_blanking_method == "distance_function");
    const amrex::Real smooth_len = m_smoothing_length * dx[2];

    // Pass 1: surface geometry, roughness and volume fraction, including
    // ghost cells so that the neighbor search in pass 2 has valid data.
    amrex::ParallelFor(
        fraction, m_terrain_fraction.num_grow(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real x = prob_lo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real y = prob_lo[1] + ((j + 0.5_rt) * dx[1]);
            const amrex::Real z = prob_lo[2] + ((k + 0.5_rt) * dx[2]);

            const auto height = [=](const amrex::Real xq,
                                    const amrex::Real yq) {
                return interp::bilinear(
                    xterrain_ptr, xterrain_ptr + xterrain_size, yterrain_ptr,
                    yterrain_ptr + yterrain_size, zterrain_ptr, xq, yq);
            };

            // Height and slopes at the cell center; slopes by central
            // differences over one cell width of the interpolated surface
            const amrex::Real terrain_ht = height(x, y);
            const amrex::Real slope_x = (height(x + 0.5_rt * dx[0], y) -
                                         height(x - 0.5_rt * dx[0], y)) /
                                        dx[0];
            const amrex::Real slope_y = (height(x, y + 0.5_rt * dx[1]) -
                                         height(x, y - 0.5_rt * dx[1])) /
                                        dx[1];
            surf_arrs[nbx](i, j, k, surf_height) = terrain_ht;
            surf_arrs[nbx](i, j, k, surf_slope_x) = slope_x;
            surf_arrs[nbx](i, j, k, surf_slope_y) = slope_y;

            // Volume fraction of terrain in the cell
            amrex::Real vol_frac = 0.0_rt;
            if (use_distance_function) {
                // Smooth blanking: 1 well below the surface, 0 well above
                const amrex::Real dist = z - terrain_ht;
                vol_frac = 0.5_rt * (1.0_rt - std::tanh(dist / smooth_len));
            } else {
                // Fraction of the cell column below the terrain height
                const amrex::Real z_bottom = prob_lo[2] + (k * dx[2]);
                vol_frac = (terrain_ht - z_bottom) / dx[2];
            }
            vol_frac = amrex::min<amrex::Real>(
                amrex::max<amrex::Real>(vol_frac, 0.0_rt), 1.0_rt);
            if (vol_frac < fraction_tol) {
                vol_frac = 0.0_rt;
            } else if (vol_frac > 1.0_rt - fraction_tol) {
                vol_frac = 1.0_rt;
            }
            // Ghost cells below the domain floor stay fluid so the bottom row
            // is not flagged as a terrain surface by the neighbor search
            frac_arrs[nbx](i, j, k, 0) = (z > prob_lo[2]) ? vol_frac : 0.0_rt;
            if (has_rate) {
                rate_arrs[nbx](i, j, k, 0) =
                    frac_arrs[nbx](i, j, k, 0) * drag_rate_solid;
            }

            // Roughness
            if (xrough_size > 0) {
                z0_arrs[nbx](i, j, k, 0) = interp::bilinear(
                    xrough_ptr, xrough_ptr + xrough_size, yrough_ptr,
                    yrough_ptr + yrough_size, z0rough_ptr, x, y);
            } else {
                z0_arrs[nbx](i, j, k, 0) = uniform_z0;
            }
        });
    amrex::Gpu::streamSynchronize();

    // Pass 2: cell classification. A fluid cell becomes a surface cell if any
    // of its six face neighbors is mostly solid; this catches the side walls
    // of steep terrain and buildings that a vertical-only search misses.
    const amrex::Real solid_threshold = m_solid_threshold;
    amrex::ParallelFor(
        fraction, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const auto& frac = frac_arrs[nbx];
            const amrex::Real f = frac(i, j, k, 0);
            int cell_mask = mask_fluid;
            if (f >= 1.0_rt) {
                cell_mask = mask_solid;
            } else if (f > 0.0_rt) {
                cell_mask = mask_surface;
            } else {
                const amrex::Real max_nb = amrex::max<amrex::Real>(
                    amrex::max<amrex::Real>(
                        frac(i - 1, j, k, 0), frac(i + 1, j, k, 0)),
                    amrex::max<amrex::Real>(
                        frac(i, j - 1, k, 0), frac(i, j + 1, k, 0)),
                    amrex::max<amrex::Real>(
                        frac(i, j, k - 1, 0), frac(i, j, k + 1, 0)));
                if (max_nb >= solid_threshold) {
                    cell_mask = mask_surface;
                }
            }
            mask_arrs[nbx](i, j, k, 0) = cell_mask;
        });
    amrex::Gpu::streamSynchronize();
}

void ImmersedTerrain::post_regrid_actions()
{
    const int nlevels = m_sim.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        initialize_fields(lev, m_sim.repo().mesh().Geom(lev));
    }
}

} // namespace kynema_sgf::immersedterrain
