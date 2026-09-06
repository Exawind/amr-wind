#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/iter_tools.H"
#include "ks_test_utils/test_utils.H"
#include "src/physics/ImmersedTerrain.H"
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
} // namespace

namespace kynema_sgf_tests {

// Face factors applied to the diffusion coefficients at the terrain interface
class ImmersedInterfaceDiffusionTest : public MeshTest
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

    void setup(const std::string& mode)
    {
        {
            amrex::ParmParse pp("ImmersedTerrain");
            pp.remove("interface_diffusion");
            pp.add("interface_diffusion", mode);
        }
        write_terrain("terrain.amrwind");
        populate_parameters();
        initialize_mesh();
        sim().pde_manager().register_icns();
        sim().init_physics();
        m_terrain = std::make_unique<Terrain>(sim());
        const int nlevels = sim().repo().num_active_levels();
        for (int lev = 0; lev < nlevels; ++lev) {
            m_terrain->initialize_fields(lev, sim().repo().mesh().Geom(lev));
        }
    }

    using Terrain = kynema_sgf::immersedterrain::ImmersedTerrain;
    std::unique_ptr<Terrain> m_terrain;
    const amrex::Real m_tol{1.0e-12_rt};
};

// dx = dy = dz = 32; plateau cells i = 14..17, solid for k = 0..2,
// partial (beta = 0.125) at k = 3.

TEST_F(ImmersedInterfaceDiffusionTest, none_declares_no_fields)
{
    setup("none");
    EXPECT_FALSE(sim().repo().field_exists("terrain_diffusion_xf"));
}

TEST_F(ImmersedInterfaceDiffusionTest, block)
{
    setup("block");
    const auto& fx = sim().repo().get_field("terrain_diffusion_xf");
    const auto& fz = sim().repo().get_field("terrain_diffusion_zf");
    // x face between fluid (13,10,1) and solid (14,10,1): fully blocked
    EXPECT_NEAR(utils::field_probe(fx, 0, 14, 10, 1), 0.0_rt, m_tol);
    // z face between solid (15,10,2) and partial (15,10,3): blocked
    EXPECT_NEAR(utils::field_probe(fz, 0, 15, 10, 3), 0.0_rt, m_tol);
    // z face between partial (15,10,3) and fluid (15,10,4): 1 - 0.125
    EXPECT_NEAR(utils::field_probe(fz, 0, 15, 10, 4), 0.875_rt, m_tol);
    // fluid/fluid faces untouched
    EXPECT_NEAR(utils::field_probe(fx, 0, 5, 5, 8), 1.0_rt, m_tol);
    EXPECT_NEAR(utils::field_probe(fz, 0, 5, 5, 8), 1.0_rt, m_tol);
}

TEST_F(ImmersedInterfaceDiffusionTest, no_slip)
{
    setup("no_slip");
    const auto& fx = sim().repo().get_field("terrain_diffusion_xf");
    const auto& fz = sim().repo().get_field("terrain_diffusion_zf");
    // Lateral wall at the face: dx / (dx/2) = 2
    EXPECT_NEAR(utils::field_probe(fx, 0, 14, 10, 1), 2.0_rt, m_tol);
    // Bottom face of the partial cell (15,10,3): centre 112 m, terrain
    // 100 m, so d1 = 12 m and the factor is 32/12
    EXPECT_NEAR(
        utils::field_probe(fz, 0, 15, 10, 3), 32.0_rt / 12.0_rt, 1.0e-10_rt);
    // Face between the partial cell (below threshold) and fluid: untouched
    EXPECT_NEAR(utils::field_probe(fz, 0, 15, 10, 4), 1.0_rt, m_tol);
    // Solid/solid and fluid/fluid faces untouched
    EXPECT_NEAR(utils::field_probe(fz, 0, 15, 10, 2), 1.0_rt, m_tol);
    EXPECT_NEAR(utils::field_probe(fx, 0, 5, 5, 8), 1.0_rt, m_tol);
}

} // namespace kynema_sgf_tests
