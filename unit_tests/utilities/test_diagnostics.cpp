#include "aw_test_utils/MeshTest.H"
#include "amr-wind/utilities/diagnostics.H"

namespace amr_wind_tests {

namespace {

void init_velocity(amr_wind::Field& velocity)
{
    const auto& mesh = velocity.repo().mesh();
    const int nlevels = velocity.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();

        for (amrex::MFIter mfi(velocity(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.validbox();
            const auto& farr = velocity(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const amrex::Real xc = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real yc = problo[1] + (j + 0.5) * dx[1];
                const amrex::Real zc = problo[2] + (k + 0.5) * dx[2];

                farr(i, j, k, 0) = 1.0 - std::pow(xc, 2.0);
                farr(i, j, k, 1) = -1.0 + std::pow(zc, 2.0);
                farr(i, j, k, 2) = 5.0 * std::cos(yc);
            });
        }
    }
}

void init_mac_velocity(
    amr_wind::Field& cc,
    amr_wind::Field& umac,
    amr_wind::Field& vmac,
    amr_wind::Field& wmac)
{
    const auto& mesh = cc.repo().mesh();
    const int nlevels = cc.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();

        for (amrex::MFIter mfi(cc(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.growntilebox(1);
            const auto& uarr = umac(lev).array(mfi);
            const auto& varr = vmac(lev).array(mfi);
            const auto& warr = wmac(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const amrex::Real x = problo[0] + i * dx[0];
                const amrex::Real yc = problo[1] + (j + 0.5) * dx[1];
                const amrex::Real zc = problo[2] + (k + 0.5) * dx[2];

                uarr(i, j, k) = 1.0 - std::pow(x, 2.0);
                varr(i, j, k) = -1.0 + std::pow(zc, 2.0);
                warr(i, j, k) = -3.0 * std::cos(yc);
            });
        }
    }
}

} // namespace

class DiagnosticsTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{24, 24, 8}};
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
            pp.addarr("is_periodic", amrex::Vector<int>{{1, 1, 0}});
        }
    }
    const amrex::Vector<amrex::Real> problo{{-5.0, -5.0, -2.0}};
    const amrex::Vector<amrex::Real> probhi{{5.0, 5.0, 2.0}};
};

TEST_F(DiagnosticsTest, Max_Vel)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& velocity = repo.declare_field("velocity", 3, 0);
    init_velocity(velocity);

    auto cc_results =
        amr_wind::diagnostics::PrintMaxVelLocations(repo, "cell-centered");

    // Check max's and min's, according to profiles
    const amrex::Real tol = 1.0e-10;
    // max(u)
    EXPECT_NEAR(cc_results[0], 1.0 - std::pow(0.5 * 10.0 / 24.0, 2.0), tol);
    // min(u)
    EXPECT_NEAR(cc_results[4], 1.0 - std::pow(11.5 * 10.0 / 24.0, 2.0), tol);
    // max(v)
    EXPECT_NEAR(cc_results[8], -1.0 + std::pow(3.5 * 4.0 / 8.0, 2.0), tol);
    // min(v)
    EXPECT_NEAR(cc_results[12], -1.0 + std::pow(0.5 * 4.0 / 8.0, 2.0), tol);
    // max(w)
    EXPECT_NEAR(cc_results[16], 5.0 * std::cos(0.5 * 10.0 / 24.0), tol);

    // Check locations (abs due to symmetry)
    EXPECT_NEAR(std::abs(cc_results[1]), 0.5 * 10.0 / 24.0, tol);
    EXPECT_NEAR(std::abs(cc_results[5]), 11.5 * 10.0 / 24.0, tol);
    EXPECT_NEAR(std::abs(cc_results[11]), 3.5 * 4.0 / 8.0, tol);
    EXPECT_NEAR(std::abs(cc_results[15]), 0.5 * 4.0 / 8.0, tol);
    EXPECT_NEAR(std::abs(cc_results[18]), 0.5 * 10.0 / 24.0, tol);
}

TEST_F(DiagnosticsTest, Max_MACvel)
{
    initialize_mesh();
    auto& repo = sim().repo();
    repo.declare_face_normal_field({"u_mac", "v_mac", "w_mac"}, 1, 1, 1);
    auto& umac = repo.get_field("u_mac");
    auto& vmac = repo.get_field("v_mac");
    auto& wmac = repo.get_field("w_mac");
    auto& cc = repo.declare_field("cc", 1, 0);
    init_mac_velocity(cc, umac, vmac, wmac);

    auto fc_results =
        amr_wind::diagnostics::PrintMaxMACVelLocations(repo, "face-centered");

    // Check max's and min's, according to profiles
    const amrex::Real tol = 1.0e-10;
    // max(umac)
    EXPECT_NEAR(fc_results[0], 1.0 - std::pow(0.0 * 10.0 / 24.0, 2.0), tol);
    // min(umac)
    EXPECT_NEAR(fc_results[4], 1.0 - std::pow(12 * 10.0 / 24.0, 2.0), tol);
    // max(vmac)
    EXPECT_NEAR(fc_results[8], -1.0 + std::pow(3.5 * 4.0 / 8.0, 2.0), tol);
    // min(vmac)
    EXPECT_NEAR(fc_results[12], -1.0 + std::pow(0.5 * 4.0 / 8.0, 2.0), tol);
    // min(wmac)
    EXPECT_NEAR(fc_results[20], -3.0 * std::cos(0.5 * 10.0 / 24.0), tol);

    // Check locations
    EXPECT_NEAR(fc_results[1], 0.0 * 10.0 / 24.0, tol);
    EXPECT_NEAR(std::abs(fc_results[5]), 12 * 10.0 / 24.0, tol);
    EXPECT_NEAR(std::abs(fc_results[11]), 3.5 * 4.0 / 8.0, tol);
    EXPECT_NEAR(std::abs(fc_results[15]), 0.5 * 4.0 / 8.0, tol);
    EXPECT_NEAR(std::abs(fc_results[22]), 0.5 * 10.0 / 24.0, tol);
}

} // namespace amr_wind_tests
