#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/iter_tools.H"
#include "ks_test_utils/test_utils.H"
#include "src/physics/ImmersedTerrain.H"
#include "src/equation_systems/icns/icns.H"
#include "src/equation_systems/icns/source_terms/ImmersedDragForcing.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

#include <limits>

using namespace amrex::literals;

namespace {
// 100 m plateau for x in [449, 576], flat ground elsewhere
void write_terrain(const std::string& fname)
{
    std::ofstream os(fname);
    os << "6\n2\n";
    os << "0.0\n448.0\n449.0\n576.0\n577.0\n1024.0\n";
    os << "0.0\n1024.0\n";
    os << "0.0\n0.0\n0.0\n0.0\n100.0\n100.0\n100.0\n100.0\n0.0\n0.0\n0.0\n0."
          "0\n";
}

void set_string(
    const std::string& prefix, const char* key, const std::string& v)
{
    amrex::ParmParse pp(prefix);
    pp.remove(key);
    pp.add(key, v);
}
} // namespace

namespace kynema_sgf_tests {

class ImmersedDragForcingTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{32, 32, 16}};
            pp.addarr("n_cell", ncell);
            pp.add("blocking_factor", 2);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> probhi{{1024.0_rt, 1024.0_rt, 512.0_rt}};
            pp.addarr("prob_hi", probhi);
        }
    }

    // Mesh, terrain, uniform velocity (10, 5, 0), dt = 0.5
    void setup()
    {
        write_terrain("terrain.amrwind");
        populate_parameters();
        initialize_mesh();
        auto& pde_mgr = sim().pde_manager();
        pde_mgr.register_icns();
        sim().init_physics();
        m_terrain = std::make_unique<Terrain>(sim());
        const int nlevels = sim().repo().num_active_levels();
        for (int lev = 0; lev < nlevels; ++lev) {
            m_terrain->initialize_fields(lev, sim().repo().mesh().Geom(lev));
        }
        sim()
            .repo()
            .get_field("velocity")
            .setVal({{10.0_rt, 5.0_rt, 0.0_rt}}, 1);
        sim().time().delta_t() = m_dt;
    }

    [[nodiscard]] amrex::Real drag_rate(const amrex::Real beta) const
    {
        const amrex::Real C = beta * m_drag_coefficient / m_dz;
        return (1.0_rt - std::exp(-C * m_dt)) / m_dt;
    }

    using Terrain = kynema_sgf::immersedterrain::ImmersedTerrain;
    std::unique_ptr<Terrain> m_terrain;
    const amrex::Real m_dt{0.5_rt};
    const amrex::Real m_dz{32.0_rt};
    const amrex::Real m_drag_coefficient{10.0_rt};
    const amrex::Real m_tol{
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt};
};

TEST_F(ImmersedDragForcingTest, laminar_drag_only)
{
    set_string("turbulence", "model", "Laminar");
    setup();
    auto& src_term = sim().pde_manager().icns().fields().src_term;
    src_term.setVal(0.0_rt);
    kynema_sgf::pde::icns::ImmersedDragForcing forcing(sim());
    forcing(0, kynema_sgf::FieldState::New, src_term(0));

    // Solid cell: full-rate relaxation of (10, 5, 0) toward zero
    EXPECT_NEAR(
        utils::field_probe(src_term, 0, 15, 10, 1, 0),
        -drag_rate(1.0_rt) * 10.0_rt, m_tol);
    EXPECT_NEAR(
        utils::field_probe(src_term, 0, 15, 10, 1, 1),
        -drag_rate(1.0_rt) * 5.0_rt, m_tol);
    EXPECT_NEAR(utils::field_probe(src_term, 0, 15, 10, 1, 2), 0.0_rt, m_tol);
    // Partial cell (beta = 0.125): fraction-weighted rate, no wall model
    EXPECT_NEAR(
        utils::field_probe(src_term, 0, 15, 10, 3, 0),
        -drag_rate(0.125_rt) * 10.0_rt, m_tol);
    // Side-wall surface cell and far fluid cell: untouched in laminar flow
    EXPECT_NEAR(utils::field_probe(src_term, 0, 13, 10, 1, 0), 0.0_rt, m_tol);
    EXPECT_NEAR(utils::field_probe(src_term, 0, 13, 10, 1, 1), 0.0_rt, m_tol);
    EXPECT_NEAR(utils::field_probe(src_term, 0, 5, 5, 8, 0), 0.0_rt, m_tol);
}

TEST_F(ImmersedDragForcingTest, turbulent_wall_model_faces)
{
    set_string("turbulence", "model", "Smagorinsky");
    set_string("ImmersedDragForcing", "wall_model", "cell_offset");
    setup();
    auto& src_term = sim().pde_manager().icns().fields().src_term;
    src_term.setVal(0.0_rt);
    kynema_sgf::pde::icns::ImmersedDragForcing forcing(sim());
    forcing(0, kynema_sgf::FieldState::New, src_term(0));

    // Solid interior: drag only, wall model is skipped at beta = 1
    EXPECT_NEAR(
        utils::field_probe(src_term, 0, 15, 10, 1, 0),
        -drag_rate(1.0_rt) * 10.0_rt, m_tol);
    // Partial cell on top of the plateau: bottom face active, wall model
    // decelerates both horizontal components beyond the drag alone
    EXPECT_LT(
        utils::field_probe(src_term, 0, 15, 10, 3, 0),
        (-drag_rate(0.125_rt) * 10.0_rt) - m_tol);
    EXPECT_LT(
        utils::field_probe(src_term, 0, 15, 10, 3, 1),
        (-drag_rate(0.125_rt) * 5.0_rt) - m_tol);
    EXPECT_NEAR(utils::field_probe(src_term, 0, 15, 10, 3, 2), 0.0_rt, m_tol);
    // Fluid cell beside the west wall of the plateau: only the east face is
    // active, so the wall-normal u is untouched and the tangential v is
    // decelerated. Nothing but the six-face search reaches this cell.
    EXPECT_NEAR(utils::field_probe(src_term, 0, 13, 10, 1, 0), 0.0_rt, m_tol);
    EXPECT_LT(utils::field_probe(src_term, 0, 13, 10, 1, 1), -m_tol);
    EXPECT_NEAR(utils::field_probe(src_term, 0, 13, 10, 1, 2), 0.0_rt, m_tol);
    // Far fluid cell
    EXPECT_NEAR(utils::field_probe(src_term, 0, 5, 5, 8, 0), 0.0_rt, m_tol);
}

TEST_F(ImmersedDragForcingTest, terrain_height_matches_normal_on_flat_top)
{
    set_string("turbulence", "model", "Smagorinsky");
    set_string("ImmersedDragForcing", "wall_model", "terrain_height");
    setup();
    auto& src_term = sim().pde_manager().icns().fields().src_term;

    src_term.setVal(0.0_rt);
    {
        kynema_sgf::pde::icns::ImmersedDragForcing forcing(sim());
        forcing(0, kynema_sgf::FieldState::New, src_term(0));
    }
    amrex::Vector<amrex::Real> height_method(AMREX_SPACEDIM);
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        height_method[n] = utils::field_probe(src_term, 0, 15, 10, 3, n);
    }

    // On the flat plateau top the slopes vanish, the normal is vertical and
    // the surface-normal method must reproduce the terrain-height method:
    // d1 = 12, d2 = 44, reference cell k+1 in both.
    set_string("ImmersedDragForcing", "wall_model", "surface_normal");
    src_term.setVal(0.0_rt);
    {
        kynema_sgf::pde::icns::ImmersedDragForcing forcing(sim());
        forcing(0, kynema_sgf::FieldState::New, src_term(0));
    }
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        EXPECT_NEAR(
            utils::field_probe(src_term, 0, 15, 10, 3, n), height_method[n],
            m_tol);
    }
    // And it differs from the crude 0.5 dz / 1.5 dz offsets
    set_string("ImmersedDragForcing", "wall_model", "cell_offset");
    src_term.setVal(0.0_rt);
    {
        kynema_sgf::pde::icns::ImmersedDragForcing forcing(sim());
        forcing(0, kynema_sgf::FieldState::New, src_term(0));
    }
    EXPECT_GT(
        std::abs(
            utils::field_probe(src_term, 0, 15, 10, 3, 0) - height_method[0]),
        m_tol);
}

} // namespace kynema_sgf_tests
