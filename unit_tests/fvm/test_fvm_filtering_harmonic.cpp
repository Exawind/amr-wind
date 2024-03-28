#include "aw_test_utils/MeshTest.H"
#include "amr-wind/fvm/filter_harmonic.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "AMReX_ParmParse.H"

namespace amr_wind_tests {

class FvmOpTestFilteringHarmonic : public MeshTest
{};

namespace {

void initialize_scalar(
    const amrex::Box& bx, const amrex::Array4<amrex::Real>& scalar_arr)
{
    // grow the box by 1 so that x,y,z go out of bounds and min(max())
    // corrects it and it fills the ghosts with wall values
    amrex::ParallelFor(grow(bx, 1), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        scalar_arr(i, j, k) = (amrex::Real)3 + i + 2 * j + 3 * k;
    });
}

amrex::Real filtering_test_impl(amr_wind::Field& scalar, const int pdegree)
{
    const int ncoeff = (pdegree + 1) * (pdegree + 1) * (pdegree + 1);

    amrex::Gpu::DeviceVector<amrex::Real> coeff(ncoeff, 0.00123);

    const auto& geom = scalar.repo().mesh().Geom();

    run_algorithm(scalar, [&](const int lev, const amrex::MFIter& mfi) {
        auto scalar_arr = scalar(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_scalar(bx, scalar_arr);
    });

    auto filtered_scalar = amr_wind::fvm::filter_harmonic(scalar);

    const int nlevels = scalar.repo().num_active_levels();
    amrex::Real error_total = 0.0;
    const amrex::Real* coeff_ptr = coeff.data();

    for (int lev = 0; lev < nlevels; ++lev) {

        error_total += amrex::ReduceSum(
            (*filtered_scalar)(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& filter_arr)
                -> amrex::Real {
                amrex::Real error = 0.0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                    const amrex::Real ph111 =
                        (amrex::Real)3 + i + 2 * j + 3 * k;
                    const amrex::Real ph011 =
                        (amrex::Real)3 + (i - 1) + 2 * j + 3 * k;
                    const amrex::Real ph211 =
                        (amrex::Real)3 + (i + 1) + 2 * j + 3 * k;
                    const amrex::Real ph101 =
                        (amrex::Real)3 + i + 2 * (j - 1) + 3 * k;
                    const amrex::Real ph121 =
                        (amrex::Real)3 + i + 2 * (j + 1) + 3 * k;
                    const amrex::Real ph110 =
                        (amrex::Real)3 + i + 2 * j + 3 * (k - 1);
                    const amrex::Real ph112 =
                        (amrex::Real)3 + i + 2 * j + 3 * (k + 1);

                    const amrex::Real inv_filx =
                        0.25 / ph011 + 0.5 / ph111 + 0.25 / ph211;
                    const amrex::Real inv_fily =
                        0.25 / ph101 + 0.5 / ph111 + 0.25 / ph121;
                    const amrex::Real inv_filz =
                        0.25 / ph110 + 0.5 / ph111 + 0.25 / ph112;

                    const amrex::Real phi_fil =
                        3.0 / (inv_filx + inv_fily + inv_filz);

                    error += std::abs(filter_arr(i, j, k) - phi_fil);
                });

                return error;
            });
    }

    return error_total;
}

} // namespace

TEST_F(FvmOpTestFilteringHarmonic, filter)
{
    // This test only evaluates the expected behavior in the bulk
    // -- relying on ordinary filter unit test for bc behavior

    constexpr double tol = 1.0e-11;

    populate_parameters();
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<int> periodic{{1, 1, 1}};
        pp.addarr("is_periodic", periodic);
    }

    initialize_mesh();

    auto& repo = sim().repo();
    const int ncomp = 1;
    const int nghost = 1;
    auto& scalar = repo.declare_field("scalar", ncomp, nghost);
    const int pdegree = 1;
    auto error_total = filtering_test_impl(scalar, pdegree);

    amrex::ParallelDescriptor::ReduceRealSum(error_total);

    EXPECT_NEAR(error_total, 0.0, tol);
}

} // namespace amr_wind_tests
