#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"

namespace amr_wind_tests {

class VOFOpTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{4, 4, 4}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 4);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
            amrex::Vector<amrex::Real> probhi{{1.0, 1.0, 1.0}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
};

namespace {

void initialize_volume_fractions(
    const int dir,
    const amrex::Box& bx,
    const amrex::Array4<amrex::Real>& vof_arr)
{
    // grow the box by 1 so that x,y,z go out of bounds and min(max()) corrects
    // it and it fills the ghosts with wall values
    const int d = dir;
    amrex::ParallelFor(grow(bx, 1), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        int ii = (d != 0 ? i : 0);
        int jj = (d != 1 ? j : 0);
        int kk = (d != 2 ? k : 0);
        if (ii + jj + kk > 3) {
            vof_arr(i, j, k) = 0.0;
        }
        if (ii + jj + kk == 3) {
            vof_arr(i, j, k) = 0.5;
        }
        if (ii + jj + kk < 3) {
            vof_arr(i, j, k) = 1.0;
        }
    });
}

void init_vof(amr_wind::Field& vof, const int dir)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        auto vof_arr = vof(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_volume_fractions(dir, bx, vof_arr);
    });
}

void initialize_volume_fractions_horizontal(
    const int dir,
    const amrex::Box& bx,
    const amrex::Real vof_val,
    const amrex::Array4<amrex::Real>& vof_arr)
{
    // grow the box by 1 so that x,y,z go out of bounds and min(max()) corrects
    // it and it fills the ghosts with wall values
    const int d = dir;
    const amrex::Real vv = vof_val;
    amrex::ParallelFor(grow(bx, 1), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        int ii = (d == 0 ? i : (d == 1 ? j : k));
        if (ii > 1) {
            vof_arr(i, j, k) = 0.0;
        }
        if (ii == 1) {
            vof_arr(i, j, k) = vv;
        }
        if (ii < 1) {
            vof_arr(i, j, k) = 1.0;
        }
    });
}

void init_vof_h(amr_wind::Field& vof, const amrex::Real vof_val, const int dir)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        auto vof_arr = vof(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_volume_fractions_horizontal(dir, bx, vof_val, vof_arr);
    });
}

amrex::Real normal_vector_test_impl(amr_wind::Field& vof, const int dir)
{
    amrex::Real error_total = 0.0;
    const int d = dir;

    for (int lev = 0; lev < vof.repo().num_active_levels(); ++lev) {

        error_total += amrex::ReduceSum(
            vof(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& vof_arr)
                -> amrex::Real {
                amrex::Real error = 0.0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                    amrex::Real mx, my, mz;
                    amr_wind::multiphase::mixed_youngs_central_normal(
                        i, j, k, vof_arr, mx, my, mz);

                    int ii = (d != 0 ? i : 0);
                    int jj = (d != 1 ? j : 0);
                    int kk = (d != 2 ? k : 0);

                    // Use L1 norm, check cells where slope is known
                    if (ii + jj + kk == 3) {
                        error += amrex::Math::abs(mx - (d != 0 ? 0.5 : 0.0));
                        error += amrex::Math::abs(my - (d != 1 ? 0.5 : 0.0));
                        error += amrex::Math::abs(mz - (d != 2 ? 0.5 : 0.0));
                    }
                });

                return error;
            });
    }
    return error_total;
}

amrex::Real fit_plane_test_impl(amr_wind::Field& vof, const int dir)
{
    amrex::Real error_total = 0.0;
    const int d = dir;

    for (int lev = 0; lev < vof.repo().num_active_levels(); ++lev) {

        error_total += amrex::ReduceSum(
            vof(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& vof_arr)
                -> amrex::Real {
                amrex::Real error = 0.0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                    int ii = (d != 0 ? i : 0);
                    int jj = (d != 1 ? j : 0);
                    int kk = (d != 2 ? k : 0);
                    // Check multiphase cells
                    if (ii + jj + kk == 3) {
                        amrex::Real mx, my, mz, alpha;
                        amr_wind::multiphase::fit_plane(
                            i, j, k, vof_arr, mx, my, mz, alpha);

                        // Check slope
                        error += amrex::Math::abs(mx - (d != 0 ? 0.5 : 0.0));
                        error += amrex::Math::abs(my - (d != 1 ? 0.5 : 0.0));
                        error += amrex::Math::abs(mz - (d != 2 ? 0.5 : 0.0));
                        // Check intercept
                        error += amrex::Math::abs(alpha - 0.5);
                    }
                });

                return error;
            });
    }
    return error_total;
}

amrex::Real fit_plane_test_impl_h(
    amr_wind::Field& vof, const amrex::Real vof_val, const int dir)
{
    amrex::Real error_total = 0.0;
    const int d = dir;
    const amrex::Real vv = vof_val;

    for (int lev = 0; lev < vof.repo().num_active_levels(); ++lev) {

        error_total += amrex::ReduceSum(
            vof(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& vof_arr)
                -> amrex::Real {
                amrex::Real error = 0.0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                    int ii = (d == 0 ? i : (d == 1 ? j : k));
                    // Check multiphase cells
                    if (ii == 1) {
                        amrex::Real mx, my, mz, alpha;
                        amr_wind::multiphase::fit_plane(
                            i, j, k, vof_arr, mx, my, mz, alpha);

                        // Check slope
                        error += amrex::Math::abs(mx - (d == 0 ? 1.0 : 0.0));
                        error += amrex::Math::abs(my - (d == 1 ? 1.0 : 0.0));
                        error += amrex::Math::abs(mz - (d == 2 ? 1.0 : 0.0));
                        // Check intercept
                        error += amrex::Math::abs(alpha - vv);
                    }
                });

                return error;
            });
    }
    return error_total;
}

} // namespace

TEST_F(VOFOpTest, volume_intercept)
{
    // Initialize random number generator
    amrex::InitRandom(0);
    for (int n = 0; n < 100; ++n) {
        amrex::Real mx = amrex::Random();
        amrex::Real my = amrex::Random();
        amrex::Real mz = amrex::Random();
        amrex::Real vof = amrex::Random();
        // Scale slope values, like in fit_plane
        amrex::Real mm2 = mx + my + mz;
        mx = mx / mm2;
        my = my / mm2;
        mz = mz / mm2;
        // Limit vof values to multiphase range
        vof = amrex::max<amrex::Real>(
            1e-12, amrex::min<amrex::Real>(1.0 - 1e-12, vof));
        // Get intercept value and check for nan
        amrex::Real alpha =
            amr_wind::multiphase::volume_intercept(mx, my, mz, vof);
        EXPECT_EQ(alpha, alpha);

        // Set one of the slope components to 0, then try again
        auto idx = (int)std::floor(amrex::Random() * 3.0);
        switch (idx) {
        case 0:
            mx = 0.0;
            break;
        case 1:
            my = 0.0;
            break;
        default:
            mz = 0.0;
            break;
        }
        // Scale slope values, like in fit_plane
        mm2 = mx + my + mz;
        mx = mx / mm2;
        my = my / mm2;
        mz = mz / mm2;
        // Get intercept value and check for nan
        alpha = amr_wind::multiphase::volume_intercept(mx, my, mz, vof);
        EXPECT_EQ(alpha, alpha);
    }

    // Check specific problem case
    amrex::Real alpha =
        amr_wind::multiphase::volume_intercept(0.5, 0.5, 0.0, 0.5);
    EXPECT_EQ(alpha, alpha);
}

TEST_F(VOFOpTest, interface_normal)
{

    constexpr double tol = 1.0e-11;

    populate_parameters();
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<int> periodic{{0, 0, 0}};
        pp.addarr("is_periodic", periodic);
    }

    initialize_mesh();

    auto& repo = sim().repo();
    const int ncomp = 1;
    const int nghost = 3;
    auto& vof = repo.declare_field("vof", ncomp, nghost);

    amrex::Real error_total = 0.0;
    // constant in x
    init_vof(vof, 0);
    error_total = normal_vector_test_impl(vof, 0);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, tol);
    // constant in y
    init_vof(vof, 1);
    error_total = normal_vector_test_impl(vof, 1);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, tol);
    // constant in z
    init_vof(vof, 2);
    error_total = normal_vector_test_impl(vof, 2);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, tol);
}

TEST_F(VOFOpTest, interface_plane)
{

    constexpr double tol = 1.0e-11;

    populate_parameters();
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<int> periodic{{0, 0, 0}};
        pp.addarr("is_periodic", periodic);
    }

    initialize_mesh();

    auto& repo = sim().repo();
    const int ncomp = 1;
    const int nghost = 3;
    auto& vof = repo.declare_field("vof", ncomp, nghost);

    /* -- Diagonal plane, 2D orientation -- */
    amrex::Real error_total = 0.0;
    // constant in x
    init_vof(vof, 0);
    error_total = fit_plane_test_impl(vof, 0);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, tol);
    // constant in y
    init_vof(vof, 1);
    error_total = fit_plane_test_impl(vof, 1);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, tol);
    // constant in z
    init_vof(vof, 2);
    error_total = fit_plane_test_impl(vof, 2);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, tol);

    /* -- "Horizontal" plane, 1D orientation -- */
    amrex::InitRandom(0);
    for (int n = 0; n < 20; ++n) {
        amrex::Real vof_val = amrex::Random();
        vof_val = amrex::max<amrex::Real>(
            1e-12, amrex::min<amrex::Real>(1.0 - 1e-12, vof_val));
        // in x
        init_vof_h(vof, vof_val, 0);
        error_total = fit_plane_test_impl_h(vof, vof_val, 0);
        amrex::ParallelDescriptor::ReduceRealSum(error_total);
        EXPECT_NEAR(error_total, 0.0, tol);
        // in y
        init_vof_h(vof, vof_val, 1);
        error_total = fit_plane_test_impl_h(vof, vof_val, 1);
        amrex::ParallelDescriptor::ReduceRealSum(error_total);
        EXPECT_NEAR(error_total, 0.0, tol);
        // in z
        init_vof_h(vof, vof_val, 2);
        error_total = fit_plane_test_impl_h(vof, vof_val, 2);
        amrex::ParallelDescriptor::ReduceRealSum(error_total);
        EXPECT_NEAR(error_total, 0.0, tol);
    }
}

} // namespace amr_wind_tests
