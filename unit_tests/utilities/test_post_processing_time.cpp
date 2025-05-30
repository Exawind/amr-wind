
#include "aw_test_utils/MeshTest.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/PostProcessing.H"

namespace amr_wind_tests {

class PostProcTimeTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        {
            amrex::ParmParse pp("time");
            pp.add("fixed_dt", -0.1);
            pp.add("init_shrink", 0.1);
            pp.add("cfl", 0.45);
            pp.add("verbose", -1);
            pp.add("regrid_interval", -1);
            pp.add("checkpoint_interval", -1);
            pp.add("plot_interval", -1);
        }
        {
            amrex::ParmParse pp("incflo");
            pp.add("post_processing", (std::string) "fnorm");
        }
        {
            amrex::ParmParse pp("fnorm");
            pp.add("type", (std::string) "FieldNorms");
        }
    }
};

TEST_F(PostProcTimeTest, time_interval)
{
    populate_parameters();
    {
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", 0.3);
        pp.add("stop_time", 5.0);
        pp.add("max_step", 100);
    }
    {
        amrex::ParmParse pp("fnorm");
        pp.add("output_time_interval", 1.0);
    }
    initialize_mesh();

    auto& m_sim = sim();
    amr_wind::PostProcessManager& post_manager = m_sim.post_manager();
    auto& time = sim().time();
    post_manager.pre_init_actions();
    post_manager.post_init_actions();

    int out_counter = 0;
    amrex::Real out_time_sum = 0.0;
    int out_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45 / 0.3, 0.0, 0.0);
        time.advance_time();
        post_manager.post_advance_work();
    }

    // Read file output for time steps, times
    std::string fname = "post_processing/fnorm00000.txt";
    std::ifstream ifh(fname, std::ios::in);
    if (!ifh.good()) {
        amrex::Abort("Cannot find file: " + fname);
    }
    int data_tstep;
    amrex::Real data_time;
    ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    while (ifh >> data_tstep) {
        ifh >> data_time;
        ++out_counter;
        out_step_sum += data_tstep;
        out_time_sum += data_time;
    }

    EXPECT_EQ(out_counter, 1 + 5);
    EXPECT_NEAR(out_time_sum, 0. + 1.2 + 2.1 + 3.0 + 4.2 + 5.1, 1e-8);
    EXPECT_EQ(out_step_sum, 4 + 7 + 10 + 14 + 17);

    // Remove file
    if (ifh.good()) {
        remove(fname.c_str());
    }
}

TEST_F(PostProcTimeTest, enforce_time_interval)
{
    populate_parameters();
    {
        amrex::ParmParse pp("time");
        pp.add("stop_time", 1.0);
        pp.add("max_step", 10);
    }
    {
        amrex::ParmParse pp("fnorm");
        pp.add("output_time_interval", 0.5);
        pp.add("enforce_output_time_dt", true);
    }
    initialize_mesh();

    auto& m_sim = sim();
    amr_wind::PostProcessManager& post_manager = m_sim.post_manager();
    auto& time = sim().time();
    post_manager.pre_init_actions();
    // Give a real value for dt so initial output happens
    time.delta_t() = 1.0;
    post_manager.post_init_actions();

    int out_counter = 0;
    amrex::Real out_time_sum = 0.0;
    int out_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45 / 0.4, 0.0, 0.0);
        time.advance_time();
        post_manager.post_advance_work();
    }

    // Read file output for time steps, times
    std::string fname = "post_processing/fnorm00000.txt";
    std::ifstream ifh(fname, std::ios::in);
    if (!ifh.good()) {
        amrex::Abort("Cannot find file: " + fname);
    }
    int data_tstep;
    amrex::Real data_time;
    ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    while (ifh >> data_tstep) {
        ifh >> data_time;
        ++out_counter;
        out_step_sum += data_tstep;
        out_time_sum += data_time;
    }

    EXPECT_EQ(out_counter, 1 + 2);
    EXPECT_NEAR(out_time_sum, 1.5, 1e-8);
    EXPECT_EQ(out_step_sum, 2 + 6);

    // Remove file
    if (ifh.good()) {
        remove(fname.c_str());
    }
}

} // namespace amr_wind_tests