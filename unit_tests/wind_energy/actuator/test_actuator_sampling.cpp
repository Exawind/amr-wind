#include "aw_test_utils/MeshTest.H"

#include "amr-wind/wind_energy/actuator/ActuatorContainer.H"
#include "amr-wind/core/gpu_utils.H"
#include "amr-wind/core/vs/vector_space.H"

#include <algorithm>

namespace amr_wind_tests {
namespace {

// Utility function to populate the velocity field used for tests
void init_field(amr_wind::Field& fld)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();
    const int ncomp = fld.num_comp();

    amrex::Real offset = 0.0;
    if (fld.field_location() == amr_wind::FieldLoc::CELL) offset = 0.5;

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();

        for (amrex::MFIter mfi(fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.growntilebox();
            const auto& farr = fld(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const amrex::Real x = problo[0] + (i + offset) * dx[0];
                const amrex::Real y = problo[1] + (j + offset) * dx[1];
                const amrex::Real z = problo[2] + (k + offset) * dx[2];

                for (int d = 0; d < ncomp; d++) farr(i, j, k, d) = x + y + z;
            });
        }
    }
}

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

class ActuatorTest : public MeshTest
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

} // namespace

TEST_F(ActuatorTest, act_container)
{
    const int nprocs = amrex::ParallelDescriptor::NProcs();
    if (nprocs > 2) {
        GTEST_SKIP();
        return;
    }

    const int iproc = amrex::ParallelDescriptor::MyProc();
    initialize_mesh();
    auto& vel = sim().repo().declare_field("velocity", 3, 3);
    init_field(vel);

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

    ac.sample_velocities(vel);

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

} // namespace amr_wind_tests
