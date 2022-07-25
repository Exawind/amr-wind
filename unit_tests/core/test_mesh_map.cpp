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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real eval_coord(
    const amrex::Real x,
    const amrex::Real beta,
    const amrex::Real prob_lo,
    const amrex::Real len)
{
    return (beta == 0.0)
               ? x
               : (prob_lo +
                  len / 2 *
                      (1 - std::tanh(beta * (1 - 2 * (x - prob_lo) / len)) /
                               std::tanh(beta)));
}

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
        Vector str_coord{{0.0, 0.0, 0.0}};

        auto uni_coord = sim().mesh_mapping()->stretched_to_unstretched(
            str_coord, m_mesh->Geom(0));

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            EXPECT_NEAR(str_coord[i], uni_coord[i], 1e-12) << i;
        }
    }
    // top unchanged
    {
        Vector str_coord{{32.0, 32.0, 32.0}};

        auto uni_coord = sim().mesh_mapping()->stretched_to_unstretched(
            str_coord, m_mesh->Geom(0));

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            EXPECT_NEAR(str_coord[i], uni_coord[i], 1e-12) << i;
        }
    }
    {
        Vector gold{{3.0, 0.5, 0.5}};

        amrex::Real x = eval_coord(gold[0], 0.0, 0.0, 32.0);
        amrex::Real y = eval_coord(gold[1], 3.0, 0.0, 32.0);
        amrex::Real z = eval_coord(gold[2], 0.0, 0.0, 32.0);

        Vector str_coord{{x, y, z}};

        auto uni_coord = sim().mesh_mapping()->stretched_to_unstretched(
            str_coord, m_mesh->Geom(0));

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {

            EXPECT_NEAR(gold[i], uni_coord[i], 1e-12) << i;
        }
    }
}

} // namespace amr_wind_tests
