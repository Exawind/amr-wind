#include "mms_test_utils.H"
#include "trig_ops.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"

#include "masa.h"
#include "MMS.H"

namespace amr_wind_tests {
TEST_F(MMSMeshTest, mms_error)
{
    populate_parameters();
    utils::populate_mms_params();

    initialize_mesh();
    auto& frepo = mesh().field_repo();
    auto& velocityf = frepo.declare_field("velocity", 3, 0);
    frepo.declare_field("density");

    amr_wind::mms::MMS mms(sim());
    const int nlevels = mesh().num_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        mms.initialize_fields(lev, mesh().Geom(lev));
    }

    // perturb the velocity fields from the initial condition
    auto velocity = velocityf.vec_ptrs();
    run_algorithm(
        mesh().num_levels(), velocity,
        [&](const int lev, const amrex::MFIter& mfi) {
            auto vel = velocity[lev]->array(mfi);

            const auto& bx = mfi.validbox();
            const auto& dx = mesh().Geom(lev).CellSizeArray();
            const auto& problo = mesh().Geom(lev).ProbLoArray();

            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                    const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                    vel(i, j, k, 0) += std::sin(amr_wind::utils::two_pi() * x);
                    vel(i, j, k, 1) +=
                        std::sin(amr_wind::utils::two_pi() * 2.0 * y);
                    vel(i, j, k, 2) += std::cos(amr_wind::utils::two_pi() * z);
                });
        });

    const amrex::Real u_mms_err =
        mms.compute_error(0, velocityf, masa_eval_3d_exact_u);
    const amrex::Real v_mms_err =
        mms.compute_error(1, velocityf, masa_eval_3d_exact_v);
    const amrex::Real w_mms_err =
        mms.compute_error(2, velocityf, masa_eval_3d_exact_w);

    const amrex::Real tol = 1.0e-12;

    amrex::Real golds[AMREX_SPACEDIM] = {
        0.67158412372224829, 0.69787170811260568, 0.74092831283749405};
    EXPECT_NEAR(u_mms_err, golds[0], tol);
    EXPECT_NEAR(v_mms_err, golds[1], tol);
    EXPECT_NEAR(w_mms_err, golds[2], tol);
}

} // namespace amr_wind_tests
