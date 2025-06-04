
#include "aw_test_utils/MeshTest.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/PostProcessing.H"
#include "aw_test_utils/test_utils.H"

namespace amr_wind_tests {

class TimeAveragingTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        {
            amrex::ParmParse pp("time");
            pp.add("fixed_dt", 0.1);
            pp.add("cfl", 0.45);
            pp.add("verbose", -1);
            pp.add("regrid_interval", -1);
            pp.add("checkpoint_interval", -1);
            pp.add("plot_interval", -1);
        }
        {
            amrex::ParmParse pp("incflo");
            pp.addarr("post_processing", amrex::Vector<std::string>{"tavg1","tavg2"});
        }
        {
            amrex::ParmParse pp("tavg1");
            pp.add("type", (std::string) "TimeAveraging");
            pp.add("labels", (std::string) "means");
            pp.add("averaging_window", (amrex::Real) 2.0);
            pp.add("averaging_interval", (amrex::Real) 2);
        }
        {
            amrex::ParmParse pp("tavg1.means");
            pp.add("fields", (std::string) "temperature");
            pp.add("averaging_type", (std::string) "ReAveraging");
        }
        {
            amrex::ParmParse pp("tavg2");
            pp.add("type", (std::string) "TimeAveraging");
            pp.add("labels", (std::string) "means");
            pp.add("averaging_window", (amrex::Real) 2.0);
            pp.add("averaging_time_interval", (amrex::Real) 0.2);
        }
        {
            amrex::ParmParse pp("tavg2.means");
            pp.add("fields", (std::string) "temperature");
            pp.add("averaging_type", (std::string) "ReAveraging");
        }
    }

protected:
    std::string m_name1 = "temperature_mean_tavg1";
    std::string m_name2 = "temperature_mean_tavg2";
};

TEST_F(TimeAveragingTest, step_vs_time)
{
    populate_parameters();
    initialize_mesh();

    auto& m_sim = sim();
    amr_wind::PostProcessManager& post_manager = m_sim.post_manager();
    auto& time = sim().time();
    auto& temp = sim().repo().declare_cc_field("temperature",1,1,1);
    temp.setVal(0.);
    post_manager.pre_init_actions();
    post_manager.post_init_actions();

    const auto& f_avg1 = sim().repo().get_field(m_name1);
    const auto& f_avg2 = sim().repo().get_field(m_name2);
    constexpr amrex::Real tol = 1e-8;

    while (time.new_timestep()) {
        // Give field a linear profile with time
        temp.setVal(10. * time.current_time());
        time.set_current_cfl(0.45 / 0.3, 0.0, 0.0);
        time.advance_time();
        post_manager.post_advance_work();
        // Get field values
        const auto max_f1 = utils::field_max(f_avg1);
        const auto max_f2 = utils::field_max(f_avg2);
        const auto min_f1 = utils::field_min(f_avg1);
        const auto min_f2 = utils::field_min(f_avg2);
        // Confirm mean fields are uniform
        EXPECT_NEAR(max_f1,min_f1,tol);
        EXPECT_NEAR(max_f2,min_f2,tol);
        // Confirm mean fields match with different methods
        EXPECT_NEAR(max_f1,max_f2,tol);
        EXPECT_NEAR(min_f1,min_f2,tol);
        // Confirm value is as expected
    }

    
}

} // namespace amr_wind_tests