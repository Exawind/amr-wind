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
        std::cout << time.new_time() << std::endl;
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

TEST_F(SimTimeTest, plt_chk_timeinterval_loop)
{
    build_simtime_params();
    {
        amrex::ParmParse pp("time");
        pp.add("regrid_interval", -1);
        pp.add("checkpoint_interval", -1);
        pp.add("plot_interval", -1);
        pp.add("checkpoint_time_interval", 4.0);
        pp.add("plot_time_interval", 1.0);
        pp.add("fixed_dt", 0.3);
        pp.add("stop_time", 5.0);
        pp.add("max_step", 100);
    }
    amr_wind::SimTime time;
    time.parse_parameters();

    int counter = 0;
    int plot_counter = 0;
    int plot_step_sum = 0;
    int chkpt_counter = 0;
    int chkpt_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45 / 0.3, 0.0, 0.0);
        ++counter;

        if (time.write_plot_file()) {
            ++plot_counter;
            plot_step_sum += counter;
        }
        if (time.write_checkpoint()) {
            ++chkpt_counter;
            chkpt_step_sum += counter;
        }
    }
    EXPECT_EQ(plot_counter, 5);
    EXPECT_EQ(plot_step_sum, 4 + 7 + 10 + 14 + 17);
    EXPECT_EQ(chkpt_counter, 1);
    EXPECT_EQ(chkpt_step_sum, 14);
}

TEST_F(SimTimeTest, enforce_dt_out)
{
    // Should not change if already correct
    amrex::Real result = get_enforced_dt_for_output(0.1, 3.9, 2.0, 1e-3);
    EXPECT_NEAR(result, 0.1, 1e-12);

    // Should not change if short of interval
    result = get_enforced_dt_for_output(0.1, 3.9 - 2.0 * 5e-4, 2.0, 1e-3);
    EXPECT_NEAR(result, 0.1, 1e-12);

    // Should not change if starting near interval
    result = get_enforced_dt_for_output(0.1, 4.0, 2.0, 1e-3);
    EXPECT_NEAR(result, 0.1, 1e-12);
    result = get_enforced_dt_for_output(0.1, 4.0 - 2.0 * 0.99e-3, 2.0, 1e-3);
    EXPECT_NEAR(result, 0.1, 1e-12);
    // Past the tolerance, will change
    result = get_enforced_dt_for_output(0.1, 4.0 - 2.0 * 1.01e-3, 2.0, 1e-3);
    EXPECT_LT(result, 0.1);

    // Shortens dt if overlapping with intervals
    result = get_enforced_dt_for_output(0.1, 3.95, 2.0, 1e-3);
    EXPECT_NEAR(result, 0.05, 1e-12);
}

TEST_F(SimTimeTest, enforce_timeinterval)
{
    build_simtime_params();
    {
        amrex::ParmParse pp("time");
        pp.add("regrid_interval", -1);
        pp.add("checkpoint_interval", -1);
        pp.add("plot_interval", -1);

        pp.add("plot_time_interval", 0.5);
        pp.add("enforce_plot_time_dt", true);

        // Default values for tolerances
        pp.add("stop_time", 1.0);
        pp.add("max_step", 10);
    }
    amr_wind::SimTime time;
    time.parse_parameters();

    int counter = 0;
    int plot_counter = 0;
    amrex::Real plot_time_sum = 0.0;
    int plot_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45 / 0.4, 0.0, 0.0);
        ++counter;
        if (time.write_plot_file()) {
            ++plot_counter;
            plot_time_sum += time.new_time();
            plot_step_sum += counter;
        }
    }
    EXPECT_EQ(plot_counter, 2);
    EXPECT_NEAR(plot_time_sum, 1.5, 1e-8);
    EXPECT_EQ(plot_step_sum, 2 + 6);
}

TEST_F(SimTimeTest, enforce_timeinterval_bigtimetol)
{
    build_simtime_params();
    {
        amrex::ParmParse pp("time");
        pp.add("regrid_interval", -1);
        pp.add("checkpoint_interval", -1);
        pp.add("plot_interval", -1);

        pp.add("plot_time_interval", 0.5);
        pp.add("enforce_plot_time_dt", true);

        // Choose values that could give issues
        pp.add("enforce_plot_dt_reltol", 1e-8);
        pp.add("plot_time_interval_reltol", 1e0);

        pp.add("stop_time", 1.0);
        pp.add("max_step", 10);
    }
    amr_wind::SimTime time;
    time.parse_parameters();

    int counter = 0;
    int plot_counter = 0;
    amrex::Real plot_time_sum = 0.0;
    int plot_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45 / 0.4, 0.0, 0.0);
        ++counter;
        if (time.write_plot_file()) {
            ++plot_counter;
            plot_time_sum += time.new_time();
            plot_step_sum += counter;
        }
    }
    EXPECT_EQ(plot_counter, 2);
    EXPECT_NEAR(plot_time_sum, 1.5, 1e-8);
    EXPECT_EQ(plot_step_sum, 2 + 6);
}

TEST_F(SimTimeTest, enforce_timeinterval_bigdttol)
{
    build_simtime_params();
    {
        amrex::ParmParse pp("time");
        pp.add("regrid_interval", -1);
        pp.add("checkpoint_interval", -1);
        pp.add("plot_interval", -1);

        pp.add("plot_time_interval", 0.5);
        pp.add("enforce_plot_time_dt", true);

        // Weak enforcement of plot time interval on dt
        pp.add("enforce_plot_dt_reltol", 1e0);
        pp.add("plot_time_interval_reltol", 1e-8);

        pp.add("stop_time", 1.0);
        pp.add("max_step", 10);
    }
    amr_wind::SimTime time;
    time.parse_parameters();

    int counter = 0;
    int plot_counter = 0;
    amrex::Real plot_time_sum = 0.0;
    int plot_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45 / 0.4, 0.0, 0.0);
        ++counter;
        if (time.write_plot_file()) {
            ++plot_counter;
            plot_time_sum += time.new_time();
            plot_step_sum += counter;
        }
    }
    // With big dt tolerance, dt never gets shortened except at the end
    // (reaching t = 1.0). Ordinary plot time interval tolerance ensures that
    // plot files are still written at the first step after interval is passed.
    EXPECT_EQ(plot_counter, 2);
    EXPECT_NEAR(plot_time_sum, 0.8 + 1.0, 1e-8);
    EXPECT_EQ(plot_step_sum, 2 + 3);
}

} // namespace amr_wind_tests
