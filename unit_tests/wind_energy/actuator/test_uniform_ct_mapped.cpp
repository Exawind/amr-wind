#include "aw_test_utils/MeshTest.H"
#include "test_act_utils.H"

#include "amr-wind/wind_energy/actuator/Actuator.H"
#include "amr-wind/wind_energy/actuator/ActuatorContainer.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"
#include "amr-wind/wind_energy/actuator/disk/uniform_ct_ops.H"

namespace amr_wind_tests {
namespace {
namespace act = amr_wind::actuator;
namespace vs = amr_wind::vs;

struct UniformCt : public act::UniformCt
{
    static std::string identifier() { return "TestUniformCt"; }
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

    void declare_fields()
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
        amrex::ParmParse pp("Actuator.TestUniformCtDisk");
        pp.add("num_points_t", 3);
        pp.add("num_points_r", 3);
        pp.add("epsilon", 1.0);
        pp.add("density", 1.0);
        pp.add("rotor_diameter", 10.0);
        pp.addarr("disk_center", amrex::Vector<amrex::Real>{16.0, 16.0, 16.0});
        pp.addarr("disk_normal", amrex::Vector<amrex::Real>{-1.0, 0.0, 0.0});
        pp.addarr("thrust_coeff", amrex::Vector<amrex::Real>{1});
    }
};
} // namespace
} // namespace amr_wind_tests

namespace amr_wind {
namespace actuator {
namespace ops {
template <>
struct ReadInputsOp<::amr_wind_tests::UniformCt, ActSrcDisk>
{
    void operator()(
        ::amr_wind_tests::UniformCt::DataType& data, const utils::ActParser& pp)
    {
        ReadInputsOp<::amr_wind::actuator::UniformCt, ActSrcDisk> actual_op;
        EXPECT_NO_FATAL_FAILURE(actual_op(data, pp));

        const auto& meta = data.meta();

        // check all the necessary arrays are populated
        EXPECT_EQ(meta.num_vel_pts / 2, meta.num_force_pts);
    }
};
} // namespace ops
} // namespace actuator
} // namespace amr_wind