#include "aw_test_utils/MeshTest.H"
#include "test_act_utils.H"

#include "amr-wind/wind_energy/actuator/Actuator.H"
#include "amr-wind/wind_energy/actuator/ActuatorContainer.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"
#include "amr-wind/wind_energy/actuator/disk/Joukowsky_ops.H"

namespace amr_wind_tests {
namespace {
namespace act = amr_wind::actuator;
namespace vs = amr_wind::vs;

struct Joukowsky : public act::Joukowsky
{
    static std::string identifier() { return "TestJoukowsky"; }
};

class ActJoukowskyTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{32, 32, 32}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 16);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
            amrex::Vector<amrex::Real> probhi{{32.0, 32.0, 32.0}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }

    void intialize_domain()
    {
        initialize_mesh();
        sim().repo().declare_field("actuator_src_term", 3, 0);
        auto& vel = sim().repo().declare_field("velocity", 3, 3);
        auto& density = sim().repo().declare_field("density", 1, 3);
        vel.setVal(10.0, 0, 1, 3);
        density.setVal(1.0);
        amr_wind::actuator::ActuatorContainer::ParticleType::NextID(1U);
    }

    static void
    add_actuators(const std::string& type, amrex::Vector<std::string> labels)
    {
        amrex::ParmParse pp("Actuator");
        pp.add("type", type);
        pp.addarr("labels", labels);
    }

    static void basic_disk_setup()
    {
        amrex::ParmParse pp("Actuator.TestJoukowskyDisk");
        pp.add("num_blades", 3);
        pp.add("num_points_t", 3);
        pp.add("num_points_r", 3);
        pp.add("epsilon", 1.0);
        pp.add("density", 1.0);
        pp.add("rotor_diameter", 10.0);
        pp.addarr("disk_center", amrex::Vector<amrex::Real>{16.0, 16.0, 16.0});
        pp.addarr("disk_normal", amrex::Vector<amrex::Real>{-1.0, 0.0, 0.0});
        pp.addarr("thrust_coeff", amrex::Vector<amrex::Real>{1});
        pp.addarr("rpm", amrex::Vector<amrex::Real>{0.0});
    }
};

} // namespace

} // namespace amr_wind_tests

namespace amr_wind::actuator {
namespace ops {

template <>
struct ReadInputsOp<::amr_wind_tests::Joukowsky, ActSrcDisk>
{
    void operator()(
        ::amr_wind_tests::Joukowsky::DataType& data, const utils::ActParser& pp)
    {
        ReadInputsOp<::amr_wind::actuator::Joukowsky, ActSrcDisk> actual_op;
        EXPECT_NO_FATAL_FAILURE(actual_op(data, pp));

        const auto& meta = data.meta();

        // check all the necessary arrays are populated
        EXPECT_FALSE(meta.angular_velocity.empty());
        EXPECT_EQ(meta.num_vel_pts / 2, meta.num_force_pts);
    }
};

template <>
struct InitDataOp<::amr_wind_tests::Joukowsky, ActSrcDisk>
{
    void operator()(::amr_wind_tests::Joukowsky::DataType& data)
    {
        InitDataOp<::amr_wind::actuator::Joukowsky, ActSrcDisk> actual_op;
        EXPECT_NO_FATAL_FAILURE(actual_op(data));
        const auto& meta = data.meta();
        const auto& grid = data.grid();

        // check all the necessary arrays are sized correctly
        EXPECT_EQ(meta.root_correction.size(), meta.num_vel_pts_r);
        EXPECT_EQ(meta.tip_correction.size(), meta.num_vel_pts_r);
        EXPECT_EQ(grid.vel.size() / 2, grid.force.size());

        for (int i = 0; i < meta.num_force_pts; ++i) {
            for (int j = 0; j < 3; j++) {
                EXPECT_DOUBLE_EQ(
                    grid.pos[i][j], grid.vel_pos[i + meta.num_force_pts][j]);
            }
        }

        for (int i = 0; i < 3; ++i) {
            ASSERT_FALSE(std::isnan(meta.coplanar_vec[i]));
        }
        ASSERT_DOUBLE_EQ(meta.coplanar_vec[0], 0.0);
        ASSERT_DOUBLE_EQ(meta.coplanar_vec[1], -1.0);
    }
};

template <>
struct ComputeForceOp<::amr_wind_tests::Joukowsky, ActSrcDisk>
{
    void operator()(::amr_wind_tests::Joukowsky::DataType& data)
    {
        const auto& meta = data.meta();
        const auto& grid = data.grid();
        for (int i = 0; i < meta.num_vel_pts; ++i) {
            EXPECT_DOUBLE_EQ(10.0, grid.vel[i].x()) << ", " << i;
            EXPECT_DOUBLE_EQ(0.0, grid.vel[i].y()) << ", " << i;
            EXPECT_DOUBLE_EQ(0.0, grid.vel[i].z()) << ", " << i;
            EXPECT_DOUBLE_EQ(1.0, grid.density[i]) << ", " << i;
        }
        ComputeForceOp<::amr_wind::actuator::Joukowsky, ActSrcDisk> actual_op;
        EXPECT_NO_FATAL_FAILURE(actual_op(data));
        for (int i = 0; i < meta.num_force_pts; ++i) {
            for (int j = 0; j < 3; ++j) {
                EXPECT_FALSE(std::isnan(grid.force[i][j]))
                    << "i,j: " << i << ", " << j;
            }
            // only x direction is guaranteed to be negative
            // y and z forces will vary based on azimuthal angle
            EXPECT_GE(0.0, grid.force[i][0]) << "i: " << i;
            EXPECT_DOUBLE_EQ(1.0, grid.density[i]);
        }
    }
};

template <>
struct ProcessOutputsOp<::amr_wind_tests::Joukowsky, ActSrcDisk>
{
    ProcessOutputsOp<::amr_wind_tests::Joukowsky, ActSrcDisk>(
        ::amr_wind_tests::Joukowsky::DataType& /**/)
    {} // cppcheck-suppress missingReturn
    void operator()(::amr_wind_tests::Joukowsky::DataType& /*data*/) {}
    void read_io_options(const utils::ActParser& /**/) {}
    void prepare_outputs(const std::string& /**/) {}
    void write_outputs(){};
};

} // namespace ops
template class ::amr_wind::actuator::
    ActModel<::amr_wind_tests::Joukowsky, ::amr_wind::actuator::ActSrcDisk>;
} // namespace amr_wind::actuator

namespace amr_wind_tests {
class ActPhysicsTest : public ::amr_wind::actuator::Actuator
{
public:
    explicit ActPhysicsTest(::amr_wind::CFDSim& sim)
        : ::amr_wind::actuator::Actuator(sim)
    {}

protected:
    void prepare_outputs() final {}
};

TEST_F(ActJoukowskyTest, parsing_operations)
{
    intialize_domain();
    basic_disk_setup();
    add_actuators("TestJoukowskyDisk", {"D1"});
    ActPhysicsTest act(sim());
    act.pre_init_actions();
}

TEST_F(ActJoukowskyTest, execution)
{
    intialize_domain();
    basic_disk_setup();
    add_actuators("TestJoukowskyDisk", {"D1"});
    ActPhysicsTest act(sim());
    act.pre_init_actions();
    act.post_init_actions();
}
} // namespace amr_wind_tests