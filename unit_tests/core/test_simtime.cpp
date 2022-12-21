/** \file test_simtime.cpp
 *
 *  Unit tests for amr_wind::SimTime
 */

#include "aw_test_utils/AmrexTest.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/core/SimTime.H"

namespace amr_wind_tests {

namespace {

void build_simtime_params()
{

    amrex::ParmParse pp("time");
    pp.add("stop_time", 2.0);
    pp.add("max_step", 10);
    pp.add("fixed_dt", -0.1);
    pp.add("init_shrink", 0.1);
    pp.add("cfl", 0.45);
    pp.add("verbose", -1);
    pp.add("regrid_interval", 3);
    pp.add("plot_interval", 1);
    pp.add("checkpoint_interval", 2);
}

} // namespace

//! Create unique namespace for this test fixture
class SimTimeTest : public AmrexTest
{};

TEST_F(SimTimeTest, init)
{
    build_simtime_params();
    amr_wind::SimTime time;
    time.parse_parameters();

    constexpr double tol = 1.0e-12;
    EXPECT_EQ(time.time_index(), 0);
    EXPECT_NEAR(time.current_time(), 0.0, tol);
    EXPECT_NEAR(time.max_cfl(), 0.45, tol);

    const double cur_cfl = 0.9;
    const double dt_new = 1.0; // which comes from "2.0 * 0.45 / cur_cfl"

    // Check that the timestep size during initialization respects the shrink
    // value
    time.set_current_cfl(cur_cfl * 0.5, 0.0, 0.0);
    EXPECT_NEAR(time.deltaT(), 0.1 * dt_new, tol);
    const double first_dt = time.deltaT();

    bool stop_sim = time.new_timestep();
    EXPECT_TRUE(stop_sim);
    // Check that the timestep growth is not greater than 10% of the last
    // timestep
    time.set_current_cfl(cur_cfl * 0.5, 0.0, 0.0);
    EXPECT_NEAR(time.deltaT(), 1.1 * first_dt, tol);
}

TEST_F(SimTimeTest, time_loop)
{
    build_simtime_params();
    amr_wind::SimTime time;
    time.parse_parameters();

    int counter = 0;
    int regrid_counter = 0;
    int plot_counter = 0;
    int chkpt_counter = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(1.125, 0.0, 0.0);
        ++counter;

        if (time.write_plot_file()) {
            ++plot_counter;
        }
        if (time.write_checkpoint()) {
            ++chkpt_counter;
        }
        if (time.do_regrid()) {
            ++regrid_counter;
        }
    }
    EXPECT_EQ(counter, 5);
    EXPECT_EQ(plot_counter, 5);
    EXPECT_EQ(chkpt_counter, 2);
    EXPECT_EQ(regrid_counter, 1);

    EXPECT_TRUE(time.write_last_checkpoint());
    EXPECT_FALSE(time.write_last_plot_file());
}

TEST_F(SimTimeTest, fixed_dt_loop)
{
    build_simtime_params();
    {
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", 0.2);
    }
    amr_wind::SimTime time;
    time.parse_parameters();

    int counter = 0;
    int regrid_counter = 0;
    int plot_counter = 0;
    int chkpt_counter = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(2.0, 0.0, 0.0);
        ++counter;

        if (time.write_plot_file()) {
            ++plot_counter;
        }
        if (time.write_checkpoint()) {
            ++chkpt_counter;
        }
        if (time.do_regrid()) {
            ++regrid_counter;
        }
    }
    EXPECT_EQ(counter, 10);
    EXPECT_EQ(plot_counter, 10);
    EXPECT_EQ(chkpt_counter, 5);
    EXPECT_EQ(regrid_counter, 3);

    EXPECT_FALSE(time.write_last_checkpoint());
    EXPECT_FALSE(time.write_last_plot_file());
}

} // namespace amr_wind_tests
