#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"

namespace amr_wind_tests {

class FvmOpTest : public MeshTest
{};

namespace {

void initialize_volume_fractions(
    const amrex::Geometry& geom,
    const amrex::Box& bx,
    const int pdegree,
    amrex::Gpu::DeviceVector<amrex::Real>& cu,
    amrex::Gpu::DeviceVector<amrex::Real>& cv,
    amrex::Gpu::DeviceVector<amrex::Real>& cw,
    amrex::Array4<amrex::Real>& vel_arr)
{
    auto problo = geom.ProbLoArray();
    auto probhi = geom.ProbHiArray();
    auto dx = geom.CellSizeArray();

    const amrex::Real* cu_ptr = cu.data();
    const amrex::Real* cv_ptr = cv.data();
    const amrex::Real* cw_ptr = cw.data();

    // grow the box by 1 so that x,y,z go out of bounds and min(max()) corrects
    // it and it fills the ghosts with wall values
    amrex::ParallelFor(grow(bx, 1), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real x = amrex::min(
            amrex::max(problo[0] + (i + 0.5) * dx[0], problo[0]), probhi[0]);
        const amrex::Real y = amrex::min(
            amrex::max(problo[1] + (j + 0.5) * dx[1], problo[1]), probhi[1]);
        const amrex::Real z = amrex::min(
            amrex::max(problo[2] + (k + 0.5) * dx[2], problo[2]), probhi[2]);

        vel_arr(i, j, k, 0) =
            analytical_function::phi_eval(pdegree, cu_ptr, x, y, z);
        vel_arr(i, j, k, 1) =
            analytical_function::phi_eval(pdegree, cv_ptr, x, y, z);
        vel_arr(i, j, k, 2) =
            analytical_function::phi_eval(pdegree, cw_ptr, x, y, z);
    });
}

amrex::Real normal_vector_test_impl(amr_wind::Field& vel, const int pdegree)
{

    const int ncoeff = (pdegree + 1) * (pdegree + 1) * (pdegree + 1);

    amrex::Gpu::DeviceVector<amrex::Real> cu(ncoeff, 0.00123);
    amrex::Gpu::DeviceVector<amrex::Real> cv(ncoeff, 0.00213);
    amrex::Gpu::DeviceVector<amrex::Real> cw(ncoeff, 0.00346);

    auto& geom = vel.repo().mesh().Geom();

    run_algorithm(vel, [&](const int lev, const amrex::MFIter& mfi) {
        auto vel_arr = vel(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_velocity(geom[lev], bx, pdegree, cu, cv, cw, vel_arr);
    });

    auto str = amr_wind::fvm::strainrate(vel);

    const int nlevels = vel.repo().num_active_levels();
    amrex::Real error_total = 0.0;

    const amrex::Real* cu_ptr = cu.data();
    const amrex::Real* cv_ptr = cv.data();
    const amrex::Real* cw_ptr = cw.data();

    for (int lev = 0; lev < nlevels; ++lev) {

        const auto& problo = geom[lev].ProbLoArray();
        const auto& dx = geom[lev].CellSizeArray();

        error_total += amrex::ReduceSum(
            (*str)(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& str_arr)
                -> amrex::Real {
                amrex::Real error = 0.0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                    const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                    error += amrex::Math::abs(
                        str_arr(i, j, k) -
                        analytical_function::strainrate(
                            pdegree, cu_ptr, cv_ptr, cw_ptr, x, y, z));
                });

                return error;
            });
    }

    return error_total;
}


} // namespace

TEST_F(VOFOpTest, normal)
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
    const int ncomp = 3;
    const int nghost = 1;
    auto& vel = repo.declare_field("vel", ncomp, nghost);

    const int pdegree = 2;
    auto error_total = strainrate_test_impl(vel, pdegree);

    amrex::ParallelDescriptor::ReduceRealSum(error_total);

    EXPECT_NEAR(error_total, 0.0, tol);
}


} // namespace amr_wind_tests
