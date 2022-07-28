#include "aw_test_utils/MeshTest.H"
#include "test_act_utils.H"

#include "amr-wind/wind_energy/actuator/ActuatorContainer.H"
#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/core/gpu_utils.H"
#include "amr-wind/core/vs/vector_space.H"
#include "amr-wind/core/vs/vectorI.H"

#include <algorithm>

namespace amr_wind_tests {
namespace {

/** Mock test object to allow access to the data struct
 */
class TestActContainer : public amr_wind::actuator::ActuatorContainer
{
public:
    TestActContainer(amrex::AmrCore& mesh, const int nobj)
        : amr_wind::actuator::ActuatorContainer(mesh, nobj)
    {}

    // Accessor for the particle data holder object
    amr_wind::actuator::ActuatorCloud& get_data_obj() { return point_data(); }
};

class BaseActuatorTest : public MeshTest
{
public:
    void initalize_mesh_and_fields()
    {
        initialize_mesh();
        auto& vel = sim().repo().declare_field("velocity", 3, 2);
        auto& density = sim().repo().declare_field("density", 1, 2);
        init_field(vel);
        density.setVal(1.0);
    }

    void init_mesh_mapping()
    {
        sim().activate_mesh_map();

        if (sim().has_mesh_mapping()) {
            int n_levels = m_mesh->num_levels();
            for (int i = 0; i < n_levels; ++i) {
                sim().mesh_mapping()->create_map(i, m_mesh->Geom(i));
            }
            auto& density = sim().repo().get_field("density");
            init_field_mapped(density);
        }
    }
};

class ActuatorTest : public BaseActuatorTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{32, 32, 64}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 16);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
            amrex::Vector<amrex::Real> probhi{{128.0, 128.0, 128.0}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
};

class ActuatorMapTest : public BaseActuatorTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> n_cell({8, 8, 8});
            pp.addarr("n_cell", n_cell);
            pp.add("max_level", 0);
            pp.add("max_grid_size", 4);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
            amrex::Vector<amrex::Real> probhi{{1.0, 1.0, 1.0}};
            amrex::Vector<int> periodic{{0, 0, 0}};

            pp.addarr("is_periodic", periodic);
            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
};

} // namespace

TEST_F(ActuatorTest, act_container)
{
    const int nprocs = amrex::ParallelDescriptor::NProcs();
    if (nprocs > 2) {
        GTEST_SKIP();
    }

    const int iproc = amrex::ParallelDescriptor::MyProc();
    initalize_mesh_and_fields();
    auto& vel = sim().repo().get_field("velocity");
    auto& density = sim().repo().get_field("density");

    // Number of turbines in an MPI rank
    const int num_turbines = 2;
    // Number of points per turbine
    const int num_nodes = 16;

    TestActContainer ac(mesh(), num_turbines);
    auto& data = ac.get_data_obj();

    for (int it = 0; it < num_turbines; ++it) {
        data.num_pts[it] = num_nodes;
    }

    ac.initialize_container();

    {
        const int lev = 0;
        int idx = 0;
        const amrex::Real dz = mesh().Geom(lev).CellSize(2);
        const amrex::Real ypos = 32.0 * (iproc + 1);
        auto& pvec = data.position;
        for (int it = 0; it < num_turbines; ++it) {
            const amrex::Real xpos = 32.0 * (it + 1);
            for (int ni = 0; ni < num_nodes; ++ni) {
                const amrex::Real zpos = (ni + 0.5) * dz;

                pvec[idx].x() = xpos;
                pvec[idx].y() = ypos;
                pvec[idx].z() = zpos;
                ++idx;
            }
        }
        ASSERT_EQ(idx, ac.num_actuator_points());
    }

    ac.update_positions();

    // Check that the update positions scattered the particles to the MPI rank
    // that contains the cell enclosing the particle location
#ifndef AMREX_USE_GPU
    ASSERT_EQ(
        ac.num_actuator_points() * nprocs, ac.NumberOfParticlesAtLevel(0));
#endif
    {
        using ParIter = amr_wind::actuator::ActuatorContainer::ParIterType;
        int counter = 0;
        int total_particles = 0;
        const int lev = 0;
        for (ParIter pti(ac, lev); pti.isValid(); ++pti) {
            ++counter;
            total_particles += pti.numParticles();
        }

        if (iproc == 0) {
            ASSERT_EQ(num_turbines * num_nodes * nprocs, total_particles);
        } else {
            ASSERT_EQ(total_particles, 0);
        }
    }

    ac.sample_fields(vel, density);
    ac.Redistribute();

    // Check to make sure that the velocity sampling gathered the particles back
    // to their original MPI ranks
#ifndef AMREX_USE_GPU
    ASSERT_EQ(
        ac.num_actuator_points() * nprocs, ac.NumberOfParticlesAtLevel(0));
#endif
    {
        using ParIter = amr_wind::actuator::ActuatorContainer::ParIterType;
        int counter = 0;
        int total_particles = 0;
        const int lev = 0;
        for (ParIter pti(ac, lev); pti.isValid(); ++pti) {
            ++counter;
            total_particles += pti.numParticles();
        }

        // All particles should have been recalled back to the same patch within
        // each MPI rank
        ASSERT_EQ(counter, 1);
        // Total number of particles should be the same as what this MPI rank
        // created
        ASSERT_EQ(num_turbines * num_nodes, total_particles);
    }

    // Check the interpolated velocity field
    {
        namespace vs = amr_wind::vs;
        constexpr amrex::Real rtol = 1.0e-12;
        amrex::Real rerr = 0.0;
        const int npts = ac.num_actuator_points();
        const auto& pvec = data.position;
        const auto& vvec = data.velocity;
        for (int ip = 0; ip < npts; ++ip) {
            const auto& pos = pvec[ip];
            const auto& pvel = vvec[ip];

            const amrex::Real vval = pos.x() + pos.y() + pos.z();
            const vs::Vector vgold{vval, vval, vval};
            rerr += vs::mag_sqr(pvel - vgold);
        }
        EXPECT_NEAR(rerr, 0.0, rtol);
    }
}

TEST_F(ActuatorMapTest, containter_constant_map)
{
    const amrex::Real scale = 2.0;
    const int nprocs = amrex::ParallelDescriptor::NProcs();
    {
        amrex::ParmParse pp("geometry");
        std::string map = "ConstantMap";
        pp.add("mesh_mapping", map);
    }
    {
        amrex::ParmParse pp("ConstantMap");
        amrex::Vector<amrex::Real> scaling{{scale, scale, scale}};
        pp.addarr("scaling_factor", scaling);
    }

    initalize_mesh_and_fields();
    init_mesh_mapping();
    ASSERT_TRUE(sim().has_mesh_mapping());

    const int num_turbines = 1;
    TestActContainer ac(mesh(), num_turbines);
    auto& data = ac.get_data_obj();
    const int num_nodes = 4;
    const amrex::Real dx = 1.0 / num_nodes;
    data.num_pts[0] = num_nodes;
    ac.initialize_container();
    amr_wind::actuator::RealList golds;

    // set position of all the particles outside the integer domain
    {
        auto& pvec = data.position;
        for (int i = 0; i < num_nodes; ++i) {
            pvec[i].x() = 1.0 + (i + 0.5) * dx;
            pvec[i].y() = 1.0 + (i + 0.5) * dx;
            pvec[i].z() = 1.0 + (i + 0.5) * dx;
            golds.push_back(pvec[i].x() + pvec[i].y() + pvec[i].z());
        }
    }

    ac.update_positions(sim().mesh_mapping());
    // TODO possibly add assert in the code for this
    ASSERT_EQ(num_nodes * nprocs, ac.TotalNumberOfParticles());

    auto& vel = sim().repo().get_field("velocity");
    auto& density = sim().repo().get_field("density");
    auto& nucc = sim().repo().get_field("non_uniform_coord_cc");

    ac.sample_fields(vel, density, nucc);
    // check interpolation
    {
        auto& dvec = data.density;
        for (int i = 0; i < num_nodes; ++i) {
            EXPECT_NEAR(golds[i], dvec[i], 1.e-12) << i;
        }
    }
    ac.Redistribute();
    ASSERT_EQ(num_nodes * nprocs, ac.TotalNumberOfParticles());
}

TEST_F(ActuatorMapTest, containter_channel_flow_map)
{
    const int nprocs = amrex::ParallelDescriptor::NProcs();
    {
        amrex::ParmParse pp("geometry");
        std::string map = "ChannelFlowMap";
        pp.add("mesh_mapping", map);
    }
    {
        amrex::ParmParse pp("ChannelFlowMap");
        const amrex::Real factor = 2.0;
        amrex::Vector<amrex::Real> scaling{{factor, 0.0, 0.0}};
        pp.addarr("beta", scaling);
    }

    initalize_mesh_and_fields();
    init_mesh_mapping();
    ASSERT_TRUE(sim().has_mesh_mapping());

    const int num_turbines = 1;
    TestActContainer ac(mesh(), num_turbines);
    auto& data = ac.get_data_obj();
    const int num_nodes = 1;
    const amrex::Real dx = 1.0 / num_nodes;
    data.num_pts[0] = num_nodes;
    ac.initialize_container();
    amr_wind::actuator::RealList golds;

    // set position of all the particles as cell centers of the mapped mesh
    /*{
        const amrex::Real ypos = dx;
        auto& pvec = data.position;
        const amrex::Real xpos = dx;
        for (int ni = 0; ni < num_nodes; ++ni) {
            const amrex::Real zpos = (ni + 0.5) * dx;

            pvec[ni].x() = xpos;
            pvec[ni].y() = ypos;
            pvec[ni].z() = zpos;
            golds.push_back(pvec[ni].x() + pvec[ni].y() + pvec[ni].z());
        }
    }*/
    {
        amr_wind::vs::Vector uniform_coords{0.1, 0.5, 0.5};
        auto& pvec = data.position;
        auto stretched =
            sim().mesh_mapping()->map(uniform_coords.data(), mesh().Geom(0));
        pvec[0].x() = stretched[0];
        pvec[0].y() = stretched[1];
        pvec[0].z() = stretched[2];
        golds.push_back(pvec[0].x() + pvec[0].y() + pvec[0].z());

        auto unstretch =
            sim().mesh_mapping()->unmap(stretched.data(), mesh().Geom(0));
        for (int i = 0; i < 3; ++i) {
            ASSERT_NEAR(uniform_coords[i], unstretch[i], 1e-12);
        }
    }

    ac.update_positions(sim().mesh_mapping());
    // TODO possibly add assert in the code for this
    ASSERT_EQ(num_nodes, ac.TotalNumberOfParticles());
    auto& vel = sim().repo().get_field("velocity");
    auto& density = sim().repo().get_field("density");
    auto& nucc = sim().repo().get_field("non_uniform_coord_cc");

    ac.sample_fields(vel, density, nucc);
    // check interpolation
    auto& dvec = data.density;
    for (int i = 0; i < num_nodes; ++i) {
        EXPECT_NEAR(golds[i], dvec[i], 1.e-12) << i;
    }
}

} // namespace amr_wind_tests
