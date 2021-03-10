#include "aw_test_utils/MeshTest.H"
#include "amr-wind/fvm/curvature.H"
#include "AnalyticalFunction.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"

namespace amr_wind_tests {

class FvmOpTestCurvature : public MeshTest
{};

namespace {

void initialize_scalar(
    const amrex::Geometry& geom,
    const amrex::Box& bx,
    const int pdegree,
    amrex::Gpu::DeviceVector<amrex::Real>& c,
    amrex::Array4<amrex::Real>& scalar_arr)
{
    auto problo = geom.ProbLoArray();
    auto probhi = geom.ProbHiArray();
    auto dx = geom.CellSizeArray();

    const amrex::Real* c_ptr = c.data();

    // grow the box by 1 so that x,y,z go out of bounds and min(max()) corrects
    // it and it fills the ghosts with wall values
    amrex::ParallelFor(grow(bx, 1), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real x = amrex::min(
            amrex::max(problo[0] + (i + 0.5) * dx[0], problo[0]), probhi[0]);
        const amrex::Real y = amrex::min(
            amrex::max(problo[1] + (j + 0.5) * dx[1], problo[1]), probhi[1]);
        const amrex::Real z = amrex::min(
            amrex::max(problo[2] + (k + 0.5) * dx[2], problo[2]), probhi[2]);

        scalar_arr(i, j, k) =
            analytical_function::phi_eval(pdegree, c_ptr, x, y, z);
    });
}

amrex::Real curvature_test_impl(amr_wind::Field& scalar, const int pdegree)
{

    const int ncoeff = (pdegree + 1) * (pdegree + 1) * (pdegree + 1);

    amrex::Gpu::DeviceVector<amrex::Real> coeff(ncoeff, 0.00123);

    const auto& geom = scalar.repo().mesh().Geom();

    run_algorithm(scalar, [&](const int lev, const amrex::MFIter& mfi) {
        auto scalar_arr = scalar(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_scalar(geom[lev], bx, pdegree, coeff, scalar_arr);
    });

    auto curv_scalar = amr_wind::fvm::curvature(scalar);

    const int nlevels = scalar.repo().num_active_levels();
    amrex::Real error_total = 0.0;

    const amrex::Real* coeff_ptr = coeff.data();

    for (int lev = 0; lev < nlevels; ++lev) {

        const auto& problo = geom[lev].ProbLoArray();
        const auto& dx = geom[lev].CellSizeArray();

        error_total += amrex::ReduceSum(
            (*curv_scalar)(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& curv_arr)
                -> amrex::Real {
                amrex::Real error = 0.0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                    const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                    error += amrex::Math::abs(
                        curv_arr(i, j, k) - analytical_function::curvature(
                                                pdegree, coeff_ptr, x, y, z));
                });

                return error;
            });
    }

    return error_total;
}

} // namespace

TEST_F(FvmOpTestCurvature, curvature)
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
    const int nghost = 1;
    auto& scalar = repo.declare_field("scalar", ncomp, nghost);

    const int pdegree = 2;
    auto error_total = curvature_test_impl(scalar, pdegree);

    amrex::ParallelDescriptor::ReduceRealSum(error_total);

    EXPECT_NEAR(error_total, 0.0, tol);
}

} // namespace amr_wind_tests
