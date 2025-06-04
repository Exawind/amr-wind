#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/physics/ForestDrag.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/utilities/sampling/FieldNorms.H"

namespace {
void write_forest(const std::string& fname)
{
    std::ofstream os(fname);
    //! Write forest
    os << "1  512 256 45 200  0.2 6 0.8 \n";
    os << "1  512 512 35 200  0.2 6 0.8 \n";
    os << "1  512 612 75 200  0.2 6 0.8 \n";
    os << "2  512 762 120 200 0.2 10 0.8 \n";
}

void write_terrain(amr_wind::Field& terrain)
{
    for (int lev = 0; lev < terrain.repo().num_active_levels(); ++lev) {
        auto arrs = terrain(lev).arrays();
        for (amrex::MFIter mfi(terrain(lev)); mfi.isValid(); ++mfi) {
            auto arr = arrs[mfi.index()]; // FIX: use mfi.index() instead of mfi
            const auto bx = mfi.validbox();
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Example: flat terrain at z=10, or sloped: arr(i,j,k) = i
                    // + j;
                    arr(i, j, k) = 10.0;
                });
        }
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace

namespace amr_wind_tests {

// Testing the terrain drag reading to ensure that terrain is properly setup
class ForestTest : public MeshTest
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
    std::string m_forest_fname = "forest.amrwind";
};

TEST_F(ForestTest, forest)
{
    // Write target wind file
    write_forest(m_forest_fname);
    populate_parameters();
    initialize_mesh();
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();
    amrex::ParmParse pp("incflo");
    amrex::Vector<std::string> physics{"forestDrag"};
    pp.addarr("physics", physics);
    amr_wind::forestdrag::ForestDrag forest_drag(sim());
    const int nlevels = sim().repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = sim().repo().mesh().Geom(lev);
        forest_drag.initialize_fields(lev, geom);
    }

    constexpr amrex::Real n_forests = 3.0;
    const auto& f_id = sim().repo().get_field("forest_id");
    const amrex::Real max_id = amr_wind::field_ops::global_max_magnitude(f_id);
    EXPECT_EQ(max_id, n_forests);

    constexpr amrex::Real expected_max_drag = 0.050285714285714288;
    const auto& f_drag = sim().repo().get_field("forest_drag");
    const amrex::Real max_drag =
        amr_wind::field_ops::global_max_magnitude(f_drag);
    EXPECT_NEAR(max_drag, expected_max_drag, amr_wind::constants::TIGHT_TOL);

    constexpr amrex::Real expected_norm_drag = 0.0030635155406915832;
    const auto norm_drag =
        amr_wind::field_norms::FieldNorms::l2_norm(f_drag, 0, false);
    EXPECT_NEAR(norm_drag, expected_norm_drag, amr_wind::constants::TIGHT_TOL);
}

TEST_F(ForestTest, forest_with_terrain)
{
    write_forest(m_forest_fname);
    populate_parameters();
    initialize_mesh();

    // Add terrain field to the repo
    auto& repo = sim().repo();
    auto& terrain =
        repo.declare_field("terrain_height", 1, 0, 1, amr_wind::FieldLoc::CELL);

    write_terrain(terrain);

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();

    amrex::ParmParse pp("incflo");
    amrex::Vector<std::string> physics{"forestDrag"};
    pp.addarr("physics", physics);

    amr_wind::forestdrag::ForestDrag forest_drag(sim());
    const int nlevels = sim().repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = sim().repo().mesh().Geom(lev);
        forest_drag.initialize_fields(lev, geom);
    }

    // Now check that forest drag and id fields are set as expected
    const auto& f_id = sim().repo().get_field("forest_id");
    const amrex::Real max_id = amr_wind::field_ops::global_max_magnitude(f_id);
    EXPECT_GE(max_id, 1.0);

    const auto& f_drag = sim().repo().get_field("forest_drag");
    const amrex::Real max_drag =
        amr_wind::field_ops::global_max_magnitude(f_drag);
    EXPECT_GT(max_drag, 0.0);

    // Optionally, check that drag is only nonzero above terrain
    const auto& geom0 = sim().repo().mesh().Geom(0);
    const auto& dx = geom0.CellSizeArray();
    const auto& prob_lo = geom0.ProbLoArray();
    const auto& drag_mf = f_drag(0);
    const auto& terrain_mf = terrain(0);

    for (amrex::MFIter mfi(drag_mf); mfi.isValid(); ++mfi) {
        const auto& drag_arr = drag_mf.const_array(mfi);
        const auto& terrain_arr = terrain_mf.const_array(mfi);
        const auto& bx = mfi.validbox();
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
                    const amrex::Real terrain_z = terrain_arr(i, j, k);
                    if (drag_arr(i, j, k) != 0.0) {
                        EXPECT_GT(z, terrain_z - 1e-12)
                            << "Nonzero drag below terrain at (" << i << ","
                            << j << "," << k << ")";
                    }
                }
            }
        }
    }
}

} // namespace amr_wind_tests
