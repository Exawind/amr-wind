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

} // namespace

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

} // namespace amr_wind_tests
