#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"

namespace amr_wind_tests {

class VOFToolTest : public MeshTest
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

    amrex::Real dx = 0.25;
};

namespace {

void initialize_levelset(
    const int shape,
    const amrex::Real deltax,
    const amrex::Box& gbx,
    const amrex::Array4<amrex::Real>& lvs_arr)
{
    const int s = shape;
    const amrex::Real dx = deltax;
    amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        switch (s) {
        case 0:
            // Horizontal line
            lvs_arr(i, j, k) = 1.7 * dx;
            break;
        case 1:
            // Parabola
            lvs_arr(i, j, k) =
                1.9 * dx + 0.1 * dx * std::pow((amrex::Real)j - 0.3, 2);
            break;
        case 2:
            // Cosine profile
            lvs_arr(i, j, k) =
                2.0 * dx *
                (1.0 +
                 std::cos(((amrex::Real)i - 1.2) / amr_wind::utils::pi()));
            break;
        }
        // Subtract from local height
        lvs_arr(i, j, k) -= dx * ((amrex::Real)k + 0.5);
    });
}

/*void initialize_volume_fractions(
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
}*/

void init_lvs(
    const int dir, const amrex::Real deltax, amr_wind::Field& levelset)
{
    run_algorithm(levelset, [&](const int lev, const amrex::MFIter& mfi) {
        auto levelset_arr = levelset(lev).array(mfi);
        const auto& bx = mfi.growntilebox();
        initialize_levelset(dir, deltax, bx, levelset_arr);
    });
}

/*void init_vof(amr_wind::Field& vof, const int dir)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        auto vof_arr = vof(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_volume_fractions(dir, bx, vof_arr);
    });
}*/

amrex::Real
levelset_to_vof_test_impl(const amrex::Real deltax, amr_wind::Field& levelset)
{
    amrex::Real error_total = 0.0;
    const amrex::Real dx = deltax;

    for (int lev = 0; lev < levelset.repo().num_active_levels(); ++lev) {

        error_total += amrex::ReduceSum(
            levelset(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& levelset_arr)
                -> amrex::Real {
                amrex::Real error = 0.0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                    amrex::Real vof = amr_wind::multiphase::levelset_to_vof(
                        i, j, k, 2.0 * dx, levelset_arr);

                    // Perform checks in multiphase cells
                    if (vof > 1e-12 && vof < 1.0 - 1e-12) {
                        // Integrate to get VOF, check error
                        amrex::Real approx_vof = amrex::min(
                            1.0,
                            amrex::max(
                                0.0, (levelset_arr(i, j, k) + 0.5 * dx) / dx));
                        error += amrex::Math::abs(approx_vof - vof);
                        // std::cout << j << " " << k << " " << approx_vof << "
                        // " << vof << " " << levelset_arr(i,j,k) << std::endl;
                    }

                    // Perform checks in single-phase cells
                    if (vof <= 1e-12) {
                        // Interface should be more than half cell away,
                        // negative levelset value
                        error +=
                            amrex::max(0.0, 0.5 * dx + levelset_arr(i, j, k));
                    }
                    if (vof >= 1.0 - 1e-12) {
                        // Interface should be more than half cell away,
                        // positive levelset value
                        error +=
                            amrex::max(0.0, 0.5 * dx - levelset_arr(i, j, k));
                    }
                });

                return error;
            });
    }
    return error_total;
}
/*
amrex::Real interface_band_test_impl(amr_wind::Field& vof, const int dir)
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
                    amrex::Real mx, my, mz, alpha;

                    int ii = (d != 0 ? i : 0);
                    int jj = (d != 1 ? j : 0);
                    int kk = (d != 2 ? k : 0);
                    // Check multiphase cells
                    if (ii + jj + kk == 3) {
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
}*/

} // namespace
/*
TEST_F(VOFToolTest, interface_band)
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
        vof = std::max(1e-12, std::min(1.0 - 1e-12, vof));
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
*/
TEST_F(VOFToolTest, levelset_to_vof)
{

    constexpr double tol = 1.0e-3;

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
    auto& levelset = repo.declare_field("levelset", ncomp, nghost);
    auto& vof = repo.declare_field("vof", ncomp, nghost);

    amrex::Real error_total = 0.0;
    // profile 0: horizontal
    init_lvs(0, dx, levelset);
    error_total = levelset_to_vof_test_impl(dx, levelset);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, 1e-12);
    //  profile 1: parabola
    init_lvs(1, dx, levelset);
    error_total = levelset_to_vof_test_impl(dx, levelset);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, 0.011);
    // profile 2: cosine
    init_lvs(2, dx, levelset);
    error_total = levelset_to_vof_test_impl(dx, levelset);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, 0.016);
}

} // namespace amr_wind_tests
