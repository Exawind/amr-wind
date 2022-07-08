#include "aw_test_utils/MeshTest.H"
#include "amr-wind/mesh_mapping_models/ChannelFlowMap.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind_tests {
namespace vs = amr_wind::vs;

class MeshMapTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        {
            amrex::ParmParse amr_pp("amr");
            amrex::Vector<int> ncell{{32, 32, 32}};
            amr_pp.add("max_level", 0);
            amr_pp.add("max_grid_size", 16);
            amr_pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse geo_pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
            amrex::Vector<amrex::Real> probhi{{32.0, 32.0, 32.0}};

            geo_pp.addarr("prob_lo", problo);
            geo_pp.addarr("prob_hi", probhi);
            const std::string map_type = "ChannelFlowMap";
            geo_pp.add("mesh_mapping", map_type);
            amrex::Vector<int> periodic{0, 0, 0};
            geo_pp.addarr("is_periodic", periodic);
        }
    }
};

using Vector = amrex::Vector<amrex::Real>;

TEST_F(MeshMapTest, stretched_to_unstretched)
{
    initialize_mesh();
    sim().activate_mesh_map();

    ASSERT_TRUE(sim().has_mesh_mapping());
    int n_levels = m_mesh->num_levels();
    for (int i = 0; i < n_levels; ++i) {
        sim().mesh_mapping()->create_map(i, m_mesh->Geom(i));
    }
    // origin unchanged
    {
        Vector uni_coord{{0.0, 0.0, 0.0}};
        Vector str_coord(uni_coord);

        sim().mesh_mapping()->stretched_to_unstretched(
            uni_coord, str_coord, m_mesh->Geom(0));

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            EXPECT_NEAR(uni_coord[i], str_coord[i], 1e-12) << i;
        }
    }
    // top unchanged
    {
        Vector uni_coord{{32.0, 32.0, 32.0}};
        Vector str_coord(uni_coord);

        sim().mesh_mapping()->stretched_to_unstretched(
            uni_coord, str_coord, m_mesh->Geom(0));

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            EXPECT_NEAR(uni_coord[i], str_coord[i], 1e-12) << i;
        }
    }
}

} // namespace amr_wind_tests
