#include "aw_test_utils/MeshTest.H"
#include "test_act_utils.H"

#include "amr-wind/wind_energy/actuator/ActuatorContainer.H"
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
            auto& vel = sim().repo().declare_field("velocity", 3, 3);
            init_field_mapped(vel);
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

TEST_F(ActuatorMapTest, act_container_mesh_mapping)
{
    const int nprocs = amrex::ParallelDescriptor::NProcs();
    const int iproc = amrex::ParallelDescriptor::MyProc();
    if (nprocs > 2) {
        GTEST_SKIP();
    }

    initalize_mesh_and_fields();
    init_mesh_mapping();
    auto& vel = sim().repo().get_field("velocity");
    auto& density = sim().repo().get_field("density");

    // Number of turbines in an MPI rank
    const int num_turbines = 1;
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

            const vs::Vector vgold{pos.x(), pos.y(), pos.z()};

            vs::mag_sqr(pvel - vgold);
            for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                EXPECT_NEAR(vgold[i], pos[i], rtol) << i;
            }
        }
        EXPECT_NEAR(rerr, 0.0, rtol);
    }
}

TEST_F(ActuatorMapTest, containter_constant_map)
{
    {
        amrex::ParmParse pp("geometry");
        std::string map = "ConstantMap";
        pp.add("mesh_mapping", map);
    }
    {
        amrex::ParmParse pp("ConstantMap");
        const amrex::Real scale = 2.0;
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

    // set position of all the particles outside the integer domain
    {
        auto& pvec = data.position;
        for (int i = 0; i < num_nodes; ++i) {
            pvec[i].x() = 1.0 + (i + 0.5) * dx;
            pvec[i].y() = 1.0 + (i + 0.5) * dx;
            pvec[i].z() = 1.0 + (i + 0.5) * dx;
            amrex::Print(-1) << "x,y,z: " << pvec[i] << std::endl;
        }
    }

    ac.update_positions();
    // TODO possibly add assert in the code for this
    ASSERT_EQ(num_nodes, ac.TotalNumberOfParticles());
}

TEST_F(ActuatorMapTest, containter_channel_flow_map)
{
    {
        amrex::ParmParse pp("geometry");
        std::string map = "ChannelFlowMap";
        pp.add("mesh_mapping", map);
    }
    {
        amrex::ParmParse pp("ChannelFlowMap");
        const amrex::Real factor = 2.0;
        amrex::Vector<amrex::Real> scaling{{factor, factor, factor}};
        pp.addarr("beta", scaling);
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

    // set position of all the particles outside the integer domain
    {
        auto& pvec = data.position;
        for (int i = 0; i < num_nodes; ++i) {
            pvec[i].x() = 1.0 + (i + 0.5) * dx;
            pvec[i].y() = 1.0 + (i + 0.5) * dx;
            pvec[i].z() = 1.0 + (i + 0.5) * dx;
            amrex::Print(-1) << "x,y,z: " << pvec[i] << std::endl;
        }
    }

    ac.update_positions();
    // TODO possibly add assert in the code for this
    ASSERT_EQ(num_nodes, ac.TotalNumberOfParticles());
}

} // namespace amr_wind_tests
