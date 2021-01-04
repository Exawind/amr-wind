#include "aw_test_utils/MeshTest.H"
#include "test_act_utils.H"

#include "amr-wind/wind_energy/actuator/turbine/fast/TurbineFast.H"
#include "amr-wind/wind_energy/actuator/turbine/fast/turbine_fast_ops.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"
#include "amr-wind/wind_energy/actuator/Actuator.H"

#define AW_ENABLE_OPENFAST_UTEST 0

namespace amr_wind_tests {
namespace {

class ActTurbineFastTest : public MeshTest
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
            amrex::Vector<amrex::Real> probhi{{128.0, 128.0, 256.0}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
        {
            amrex::ParmParse pp("Actuator");
            pp.add("type", std::string("TurbineFastLine"));
        }
        {
            amrex::ParmParse pp("Actuator.TurbineFastLine");
            pp.addarr(
                "base_position", amrex::Vector<amrex::Real>{64.0, 64.0, 0.0});
            pp.add("rotor_diameter", 126.0);
            pp.add("hub_height", 90.0);
            pp.add("num_points_blade", 5);
            pp.add("num_points_tower", 5);
            pp.addarr("epsilon", amrex::Vector<amrex::Real>{3.0, 3.0, 3.0});
            pp.addarr(
                "epsilon_chord", amrex::Vector<amrex::Real>{3.0, 3.0, 3.0});

            pp.add("openfast_input_file", std::string("fast_inp/nrel5mw.fst"));
            pp.add("openfast_start_time", 0.0);
            pp.add("openfast_stop_time", 0.625);

            pp.add("nacelle_drag_coeff", 1.0);
            pp.add("nacelle_area", 100.0);
        }
    }
};

class ActTurbPhyTest : public ::amr_wind::actuator::Actuator
{
public:
    ActTurbPhyTest(::amr_wind::CFDSim& sim)
        : ::amr_wind::actuator::Actuator(sim)
    {}

protected:
    void prepare_outputs() override {}
};

} // namespace

TEST_F(ActTurbineFastTest, test_ops)
{
    namespace act = ::amr_wind::actuator;

    initialize_mesh();
    act::utils::ActParser pp("Actuator.TurbineFastLine", "Actuator.T1");
    act::ActDataHolder<act::TurbineFast> data(sim(), "T1", 0);
    {
        using ReadOp =
            act::ops::ReadInputsOp<act::TurbineFast, act::ActSrcLine>;
        ReadOp op;
        op(data, pp);
    }

    amrex::Vector<int> act_proc_count(::amrex::ParallelDescriptor::NProcs(), 0);
    act::ops::determine_root_proc<act::TurbineFast>(data, act_proc_count);

#if AW_ENABLE_OPENFAST_UTEST
    {
        using InitOp = act::ops::InitDataOp<act::TurbineFast, act::ActSrcLine>;
        InitOp op;
        op(data);
    }

    EXPECT_EQ(data.meta().num_blades, 3);
    {
        // Check all slice views
        const auto& grid = data.grid();
        const auto& blades = data.meta().blades;
        const auto& tower = data.meta().tower;
        EXPECT_EQ(std::distance(grid.pos.data(), blades[0].pos.data()), 1);
        EXPECT_EQ(std::distance(grid.pos.data(), tower.pos.data()), 16);
        EXPECT_EQ(blades[0].pos.end(), blades[1].pos.begin());
        EXPECT_EQ(blades[1].pos.end(), blades[2].pos.begin());
        EXPECT_EQ(blades[2].pos.end(), tower.pos.begin());

        EXPECT_EQ(std::distance(grid.force.data(), blades[0].force.data()), 1);
        EXPECT_EQ(std::distance(grid.force.data(), tower.force.data()), 16);
        EXPECT_EQ(blades[0].force.end(), blades[1].force.begin());
        EXPECT_EQ(blades[1].force.end(), blades[2].force.begin());
        EXPECT_EQ(blades[2].force.end(), tower.force.begin());
    }
#else
    GTEST_SKIP();
#endif
}

TEST_F(ActTurbineFastTest, fast_turbine)
{
    initialize_mesh();
    auto& vel = sim().repo().declare_field("velocity", 3, 3);
    vel.setVal(10.0, 0, 1, 3);

    {
        amrex::ParmParse pp("Actuator");
        amrex::Vector<std::string> actuators{"T1"};
        pp.addarr("labels", actuators);
    }

    ActTurbPhyTest act(sim());
    act.pre_init_actions();
#if AW_ENABLE_OPENFAST_UTEST
    act.post_init_actions();
#else
    GTEST_SKIP();
#endif
}

} // namespace amr_wind_tests
