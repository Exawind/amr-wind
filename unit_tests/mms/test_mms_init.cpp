#include "mms_test_utils.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"

#include "amr-wind/physics/mms/MMS.H"

namespace amr_wind_tests {
TEST_F(MMSMeshTest, mms_initialization)
{
#if defined(AMREX_USE_HIP)
    GTEST_SKIP();
#else
    if (!pp_utils::has_managed_memory()) GTEST_SKIP();

    populate_parameters();
    utils::populate_mms_params();

    initialize_mesh();
    auto& frepo = mesh().field_repo();
    auto& velocityf = frepo.declare_field("velocity", 3, 0);
    auto& densityf = frepo.declare_field("density");

    amr_wind::mms::MMS mms(sim());
    const int nlevels = mesh().num_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        mms.initialize_fields(lev, mesh().Geom(lev));
    }

    auto velocity = velocityf.vec_ptrs();
    auto density = densityf.vec_ptrs();
    const amrex::Real tol = 1.0e-12;

    // Test density
    {
        amrex::Real min_rho, max_rho;
        utils::field_minmax(nlevels, density, min_rho, max_rho);
        EXPECT_NEAR(min_rho, 1.0, tol);
        EXPECT_NEAR(max_rho, 1.0, tol);
    }

    // Test velocity
    {
        amrex::Vector<amrex::Real> min_vel(3), max_vel(3);
        utils::field_minmax(nlevels, velocity, min_vel, max_vel);
        EXPECT_NEAR(min_vel[0], -0.085642667425670935, tol);
        EXPECT_NEAR(min_vel[1], -0.15852969775284545, tol);
        EXPECT_NEAR(min_vel[2], -0.28208275785136977, tol);
        EXPECT_NEAR(max_vel[0], 0.085642667425670921, tol);
        EXPECT_NEAR(max_vel[1], 0.15852969775284548, tol);
        EXPECT_NEAR(max_vel[2], 0.28208275785136988, tol);
    }
#endif
}

} // namespace amr_wind_tests
