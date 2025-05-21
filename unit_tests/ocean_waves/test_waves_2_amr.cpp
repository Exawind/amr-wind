#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/ocean_waves/utils/wave_utils_K.H"
#include "amr-wind/ocean_waves/OceanWaves.H"
#include "amr-wind/ocean_waves/relaxation_zones/waves2amr_ops.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"

namespace amr_wind_tests {

class OceanWavesW2ATest : public MeshTest
{};

// This is to evaluate aspects of W2A methods that do not require wave data
TEST_F(OceanWavesW2ATest, time_advance)
{
    // Variables modified by routine
    int ntime = 0;
    bool read_flag = false;
    bool resize_flag = false;
    amrex::Real t = 0.0;
    amrex::Real t_last = -1.0;
    // Offset from wave data
    const amrex::Real t_winit = 0.0;

    // Upon initialization with waves
    int new_ntime = 1;
    amrex::Real newtime = 0.001;
    amrex::Real dt_modes = 0.1;
    t_last = 0.0; // Set within InitDataOp
    amrex::Real t_last_prior = t_last;
    int double_data = evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    // Time interpolation factor
    amrex::Real f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    // Should regrid and read, time interpolation should be ordinary. Read once
    EXPECT_EQ(double_data, 0);
    EXPECT_EQ(ntime, new_ntime);
    EXPECT_TRUE(read_flag);
    EXPECT_TRUE(resize_flag);
    EXPECT_EQ(t, 0.1);
    EXPECT_EQ(t_last, newtime);
    EXPECT_NEAR(f_interp, newtime / t, 1e-10);

    // Upon initialization without waves, when t_last has not been modified
    read_flag = false;
    resize_flag = false;
    ntime = 0;
    newtime = 0.001;
    t_last = -1.0;
    double_data = evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    EXPECT_EQ(double_data, 1);
    EXPECT_EQ(ntime, new_ntime);
    EXPECT_TRUE(read_flag);
    EXPECT_TRUE(resize_flag);
    EXPECT_EQ(t, 0.1);
    EXPECT_EQ(t_last, newtime);
    // Double data routine changes t_last_prior
    t_last_prior = (ntime - 1) * dt_modes;
    f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    EXPECT_NEAR(f_interp, newtime / t, 1e-10);

    // Upon initialization as restart, when t_last has not been modified
    read_flag = false;
    resize_flag = false;
    t_last = -1.0;
    ntime = 0;
    new_ntime = 40;
    newtime = 3.99;
    double_data = evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    EXPECT_EQ(double_data, 1);
    EXPECT_EQ(ntime, 40);
    EXPECT_TRUE(read_flag);
    EXPECT_TRUE(resize_flag);
    EXPECT_EQ(t, 40 * dt_modes);
    EXPECT_EQ(t_last, newtime);
    // Double data routine changes t_last_prior
    t_last_prior = (ntime - 1) * dt_modes;
    f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    EXPECT_NEAR(f_interp, (3.99 - 3.9) / dt_modes, 1e-10);

    // Upon initialization as restart, when newtime happens to be at ntime
    read_flag = false;
    resize_flag = false;
    t_last = -1.0;
    t_last_prior = t_last;
    ntime = 0;
    new_ntime = 40;
    newtime = 4.0;
    double_data = evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    EXPECT_EQ(double_data, 2);
    EXPECT_EQ(ntime, 40);
    EXPECT_TRUE(read_flag);
    EXPECT_TRUE(resize_flag);
    EXPECT_EQ(t, 40 * dt_modes);
    EXPECT_EQ(t_last, newtime);
    EXPECT_NEAR(f_interp, 1.0, 1e-10);

    // The following 3 tests are for during initial iterations (after the first)
    // -- Beginning from time = 0. --
    read_flag = false;
    resize_flag = false;
    ntime = 1;
    t = ntime * dt_modes;
    new_ntime = 1;
    newtime = 0.001;
    t_last = 0.001;
    t_last_prior = t_last;
    evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    // ---- Should not regrid nor read, time interpolation should do nothing
    EXPECT_EQ(ntime, new_ntime);
    EXPECT_FALSE(read_flag);
    EXPECT_FALSE(resize_flag);
    EXPECT_EQ(t_last, newtime);
    EXPECT_NEAR(f_interp, 0.0, 1e-10);
    // -- Beginning from time = 3.98 --
    read_flag = false;
    resize_flag = false;
    t_last = 3.99;
    t_last_prior = t_last;
    ntime = 40;
    t = ntime * dt_modes;
    new_ntime = 40;
    newtime = 3.99;
    evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    // ---- Should not regrid nor read, time interpolation should do nothing
    EXPECT_EQ(ntime, new_ntime);
    EXPECT_FALSE(read_flag);
    EXPECT_FALSE(resize_flag);
    EXPECT_EQ(t_last, newtime);
    EXPECT_NEAR(f_interp, 0.0, 1e-10);
    // -- Beginning from time = 3.99 --
    read_flag = false;
    resize_flag = false;
    t_last = 4.00;
    t_last_prior = t_last;
    ntime = 40;
    t = ntime * dt_modes;
    new_ntime = 40;
    newtime = 4.00;
    evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    f_interp = (newtime - t_last_prior) / (t - t_last_prior + 1e-16);
    // Should not regrid nor read, time interpolation should do nothing
    EXPECT_EQ(ntime, new_ntime);
    EXPECT_FALSE(read_flag);
    EXPECT_FALSE(resize_flag);
    EXPECT_EQ(t_last, newtime);
    EXPECT_NEAR(f_interp, 0.0, 1e-10);

    // The following tests are for during ordinary timesteps
    // -- No resize, no reading, just interp
    resize_flag = false;
    read_flag = false;
    ntime = 40;
    newtime = 3.99;
    t = 40 * dt_modes;
    t_last = 3.98;
    new_ntime = 40;
    t_last_prior = t_last;
    evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    EXPECT_EQ(ntime, 40);
    EXPECT_FALSE(read_flag);
    EXPECT_FALSE(resize_flag);
    EXPECT_EQ(t_last, newtime);
    EXPECT_NEAR(f_interp, 0.01 / 0.02, 1e-10);
    // -- No resize, yes reading
    resize_flag = false;
    read_flag = false;
    ntime = 40;
    new_ntime = 41;
    newtime = 4.01;
    t_last = 4.0;
    t_last_prior = t_last;
    evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    EXPECT_EQ(ntime, 41);
    EXPECT_TRUE(read_flag);
    EXPECT_FALSE(resize_flag);
    EXPECT_EQ(t, 41 * dt_modes);
    EXPECT_EQ(t_last, newtime);
    EXPECT_NEAR(f_interp, 0.01 / (4.1 - 4.0), 1e-10);
    // -- No resize, yes reading, big timestep means double read
    resize_flag = false;
    read_flag = false;
    ntime = 40;
    new_ntime = 42;
    newtime = 4.15;
    t_last = 4.0;
    // built-in change to t_last_prior
    t_last_prior = (new_ntime - 1) * dt_modes;
    double_data = evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    EXPECT_EQ(double_data, 1);
    EXPECT_EQ(ntime, 42);
    EXPECT_TRUE(read_flag);
    EXPECT_FALSE(resize_flag);
    EXPECT_EQ(t, 42 * dt_modes);
    EXPECT_EQ(t_last, newtime);
    // interp is from 4.1 (replaced) to 4.15
    EXPECT_NEAR(f_interp, (4.15 - 4.1) / (4.2 - 4.1), 1e-10);
    // -- No resize, yes reading, big timestep leads to convenient single read
    resize_flag = false;
    read_flag = false;
    ntime = 40;
    new_ntime = 42;
    newtime = 4.2;
    t_last = 4.0;
    t_last_prior = t_last;
    double_data = evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    EXPECT_EQ(double_data, 2);
    EXPECT_EQ(ntime, 42);
    EXPECT_TRUE(read_flag);
    EXPECT_FALSE(resize_flag);
    EXPECT_EQ(t, 42 * dt_modes);
    EXPECT_EQ(t_last, newtime);
    // interp is from 4.0 to 4.2, but it leads to 1
    EXPECT_NEAR(f_interp, (4.2 - 4.0) / (4.2 - 4.0), 1e-10);
}

TEST_F(OceanWavesW2ATest, time_reset)
{
    // Testing the arithmetic for resetting data back to the beginning
    int nstop = 10;
    int n0 = 1;

    int ntime = 11;
    int noff = update_offset_timestep(ntime, n0);
    EXPECT_EQ(ntime + noff, 1);

    // As if continued through another set of data
    ntime = nstop + ntime;
    noff = update_offset_timestep(ntime, n0);
    EXPECT_EQ(ntime + noff, 1);
}

TEST_F(OceanWavesW2ATest, time_interp)
{
    populate_parameters();
    initialize_mesh();
    // Set up fields and initial values
    auto& lvs = sim().repo().declare_cc_field("ow_levelset", 1, 3, 1);
    auto& vel = sim().repo().declare_cc_field("ow_velocity", 3, 3, 1);
    auto& w2a_lvs = sim().repo().declare_cc_field("w2a_levelset", 1, 3, 1);
    auto& w2a_vel = sim().repo().declare_cc_field("w2a_velocity", 3, 3, 1);
    lvs.setVal(0.);
    amrex::Real t_last = 0.;

    // First time from "file"
    w2a_lvs.setVal(1.0);
    amrex::Real t_w2a = 1.;

    // Interpolate once in time
    amrex::Real t_sim = 0.5;
    time_interpolate_wave_fields(
        lvs, vel, w2a_lvs, w2a_vel, t_last, t_sim, t_w2a);

    // Record time that data has been interpolated to
    t_last = t_sim;

    // Check value of levelset
    const amrex::Real fmax1 = utils::field_max(lvs);
    const amrex::Real fmin1 = utils::field_min(lvs);
    EXPECT_NEAR(fmax1, 0.5, 1e-10);
    EXPECT_NEAR(fmin1, 0.5, 1e-10);

    // Step forward in time, after current w2a data
    t_sim += 1.;

    if (t_sim > t_w2a) {
        // Copy current w2a data and update old time
        amr_wind::field_ops::copy(lvs,w2a_lvs,0,0,1,lvs.num_grow());
        t_last = t_w2a;
        // Update w2a values (as in a file read)
        w2a_lvs.setVal(0.);
        t_w2a += 1.;
    }

    // Interpolate again
    time_interpolate_wave_fields(
        lvs, vel, w2a_lvs, w2a_vel, t_last, t_sim, t_w2a);

    // Check value of levelset
    const amrex::Real fmax2 = utils::field_max(lvs);
    const amrex::Real fmin2 = utils::field_min(lvs);
    EXPECT_NEAR(fmax2, 0.5, 1e-10);
    EXPECT_NEAR(fmin2, 0.5, 1e-10);
}

} // namespace amr_wind_tests
