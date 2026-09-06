#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/iter_tools.H"
#include "ks_test_utils/test_utils.H"
#include "src/physics/ImmersedTerrain.H"
#include "src/equation_systems/icns/icns.H"
#include "src/equation_systems/temperature/source_terms/ImmersedDragTempForcing.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

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

class ImmersedDragTempForcingTest : public MeshTest
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

    // Mesh, terrain, velocity (10, 5, 0), temperature 305 K, soil 300 K
    void setup()
    {
        write_terrain("terrain.amrwind");
        populate_parameters();
        initialize_mesh();
        auto& pde_mgr = sim().pde_manager();
        pde_mgr.register_icns();
        pde_mgr.register_transport_pde("Temperature");
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
        sim().repo().get_field("temperature").setVal(m_theta, 0, 1, 1);
        sim().time().delta_t() = m_dt;
    }

    amrex::Real relax_rate(const amrex::Real beta) const
    {
        const amrex::Real C = beta * m_drag_coefficient / m_dz;
        return (1.0_rt - std::exp(-C * m_dt)) / m_dt;
    }

    amrex::Real run(const std::string& condition)
    {
        set_string("ImmersedDragTempForcing", "surface_condition", condition);
        auto& src_term =
            sim().pde_manager()("Temperature-Godunov").fields().src_term;
        src_term.setVal(0.0_rt);
        kynema_sgf::pde::temperature::ImmersedDragTempForcing forcing(sim());
        forcing(0, kynema_sgf::FieldState::New, src_term(0));
        return utils::field_probe(src_term, 0, 15, 10, 3, 0);
    }

    using Terrain = kynema_sgf::immersedterrain::ImmersedTerrain;
    std::unique_ptr<Terrain> m_terrain;
    const amrex::Real m_dt{0.5_rt};
    const amrex::Real m_dz{32.0_rt};
    const amrex::Real m_drag_coefficient{10.0_rt};
    const amrex::Real m_theta{305.0_rt};
    const amrex::Real m_soil{300.0_rt};
    const amrex::Real m_tol{1.0e-12_rt};
};

TEST_F(ImmersedDragTempForcingTest, laminar_relaxation_only)
{
    set_string("turbulence", "model", "Laminar");
    setup();
    auto& src_term =
        sim().pde_manager()("Temperature-Godunov").fields().src_term;
    src_term.setVal(0.0_rt);
    kynema_sgf::pde::temperature::ImmersedDragTempForcing forcing(sim());
    forcing(0, kynema_sgf::FieldState::New, src_term(0));

    const amrex::Real dtheta = m_theta - m_soil;
    // Solid interior relaxes toward the soil temperature at the full rate
    EXPECT_NEAR(
        utils::field_probe(src_term, 0, 15, 10, 1, 0),
        -relax_rate(1.0_rt) * dtheta, m_tol);
    // Partial cell: fraction-weighted rate, no wall model in laminar flow
    EXPECT_NEAR(
        utils::field_probe(src_term, 0, 15, 10, 3, 0),
        -relax_rate(0.125_rt) * dtheta, m_tol);
    // Side-wall surface cell and far fluid untouched
    EXPECT_NEAR(utils::field_probe(src_term, 0, 13, 10, 1, 0), 0.0_rt, m_tol);
    EXPECT_NEAR(utils::field_probe(src_term, 0, 5, 5, 8, 0), 0.0_rt, m_tol);
}

TEST_F(ImmersedDragTempForcingTest, surface_conditions)
{
    set_string("turbulence", "model", "Smagorinsky");
    set_string("ImmersedDragForcing", "wall_model", "terrain_height");
    setup();
    const amrex::Real interior = -relax_rate(0.125_rt) * (m_theta - m_soil);

    // Neutral prescribed-L: theta* = 0 with the default (huge) Obukhov length,
    // the target equals the reference temperature (uniform), so the wall
    // model adds nothing and only the interior relaxation remains
    EXPECT_NEAR(run("obukhov_length"), interior, 1.0e-9_rt);

    // Surface held at 300 K under 305 K air: heat leaves the fluid, the wall
    // model cools the cell beyond the interior relaxation
    const amrex::Real surf_temp = run("surface_temperature");
    EXPECT_LT(surf_temp, interior - 1.0e-6_rt);

    // Prescribed downward heat flux (surface cooling) must also cool
    set_string("ImmersedDragTempForcing", "surface_heat_flux", "-0.05");
    EXPECT_LT(run("heat_flux"), interior - 1.0e-6_rt);
    // and a heating flux must warm relative to the interior relaxation
    set_string("ImmersedDragTempForcing", "surface_heat_flux", "0.05");
    EXPECT_GT(run("heat_flux"), interior + 1.0e-6_rt);
}

TEST_F(ImmersedDragTempForcingTest, side_wall_cell_is_forced)
{
    set_string("turbulence", "model", "Smagorinsky");
    set_string("ImmersedDragForcing", "wall_model", "cell_offset");
    set_string(
        "ImmersedDragTempForcing", "surface_condition", "surface_temperature");
    setup();
    auto& src_term =
        sim().pde_manager()("Temperature-Godunov").fields().src_term;
    src_term.setVal(0.0_rt);
    kynema_sgf::pde::temperature::ImmersedDragTempForcing forcing(sim());
    forcing(0, kynema_sgf::FieldState::New, src_term(0));
    // Fluid cell beside the plateau wall: only the six-face search reaches it;
    // the wall at 300 K cools 305 K air
    EXPECT_LT(utils::field_probe(src_term, 0, 13, 10, 1, 0), -1.0e-6_rt);
    EXPECT_NEAR(utils::field_probe(src_term, 0, 5, 5, 8, 0), 0.0_rt, m_tol);
}

} // namespace kynema_sgf_tests
