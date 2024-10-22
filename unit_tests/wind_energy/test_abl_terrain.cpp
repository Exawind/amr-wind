#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/physics/TerrainDrag.H"

namespace {
void write_terrain(const std::string& fname)
{
    std::ofstream os(fname);
    // Write terrain height
    os << "0.0\t0.0\t0.0\n";
    os << "0.0\t1024.0\t0.0\n";
    os << "448.0\t0.0\t0.0\n";
    os << "448.0\t1024.0\t0.0\n";
    os << "449.0\t0.0\t100.0\n";
    os << "449.0\t1024.0\t100.0\n";
    os << "576.0\t0.0\t100.0\n";
    os << "576.0\t1024.0\t100.0\n";
    os << "577.0\t0.0\t0.0\n";
    os << "577.0\t1024.0\t0.0\n";
    os << "1024.0\t0.0\t0.0\n";
    os << "1024.0\t1024.0\t0.0\n";
}

} // namespace

namespace amr_wind_tests {

// Testing the terrain drag reading to ensure that terrain is properly setup
class TerrainTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        // Make computational domain like ABL mesh
        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{32, 32, 16}};
            pp.addarr("n_cell", ncell);
            pp.add("blocking_factor", 2);
        }

        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> probhi{{1024, 1024, 512}};
            pp.addarr("prob_hi", probhi);
        }
    }
    std::string terrain_fname = "terrain.amrwind";
};

TEST_F(TerrainTest, terrain)
{
    constexpr amrex::Real tol = 0;
    // Write target wind file
    write_terrain(terrain_fname);
    populate_parameters();
    initialize_mesh();
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();
    amrex::ParmParse pp("incflo");
    amrex::Vector<std::string> physics{"terrainDrag"};
    pp.addarr("physics", physics);
    amr_wind::terraindrag::TerrainDrag terrain_drag(sim());
    terrain_drag.post_init_actions();
    int value = 100;
    // Outside Point
    const auto& terrain_blank = sim().repo().get_int_field("terrain_blank");

    const int value = field_probe(terrain_blank, 0, 5, 5, 1);
    EXPECT_EQ(value, tol);
    // Inside Point
    value = terrain_drag.return_blank_value(15, 10, 1);
    EXPECT_EQ(value, 1 + tol);
}

} // namespace amr_wind_tests
