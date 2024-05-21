#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/overset/overset_ops_K.H"

namespace amr_wind_tests {

class VOFOversetOps : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{8, 8, 8}};
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
void init_vof_etc(
    amr_wind::Field& vof,
    amr_wind::Field& tg_vof,
    amr_wind::Field& norm,
    const int& dir)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        auto vof_arr = vof(lev).array(mfi);
        // tgvof and norm values are set randomly to make sure they are involved
        auto tgvof_arr = tg_vof(lev).array(mfi);
        auto norm_arr = norm(lev).array(mfi);
        const auto& bx = mfi.validbox();
        amrex::ParallelFor(
            grow(bx, 2), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                vof_arr(i, j, k) = 0.0;
                tgvof_arr(i, j, k) = 0.;
                norm_arr(i, j, k, dir) = 1.0;
                const int idx = (dir == 0 ? i : (dir == 1 ? j : k));
                if (idx == -1) {
                    // Within margin
                    vof_arr(i, j, k) = 0.41;
                    tgvof_arr(i, j, k) = 0.5;
                    norm_arr(i, j, k, dir) = -1.0;
                } else if (idx == 0) {
                    // Within margin
                    vof_arr(i, j, k) = 0.55;
                    tgvof_arr(i, j, k) = 0.5;
                    norm_arr(i, j, k, dir) = 1.0;
                } else if (idx == 1) {
                    // Outside margin, low
                    vof_arr(i, j, k) = 0.1;
                    tgvof_arr(i, j, k) = 0.2;
                    norm_arr(i, j, k, dir) = -1.0;
                } else if (idx == 2) {
                    // Also low
                    vof_arr(i, j, k) = 0.2;
                    tgvof_arr(i, j, k) = 0.22;
                    norm_arr(i, j, k, dir) = 1.0;
                } else if (idx == 3) {
                    // Above half
                    vof_arr(i, j, k) = 0.7;
                    tgvof_arr(i, j, k) = 0.8;
                    norm_arr(i, j, k, dir) = -1.0;
                } else if (idx == 4) {
                    // Also high
                    vof_arr(i, j, k) = 0.65;
                    tgvof_arr(i, j, k) = 0.91;
                    norm_arr(i, j, k, dir) = 1.0;
                } else if (idx == 5) {
                    // Also high, positive gradient
                    vof_arr(i, j, k) = 0.8;
                    tgvof_arr(i, j, k) = 0.75;
                    norm_arr(i, j, k, dir) = -1.0;
                } else if (idx == 6) {
                    // Within margin
                    vof_arr(i, j, k) = 0.45;
                    tgvof_arr(i, j, k) = 0.5;
                    norm_arr(i, j, k, dir) = 1.0;
                } else if (idx == 7) {
                    // Low, negative gradient
                    vof_arr(i, j, k) = 0.3;
                    tgvof_arr(i, j, k) = 0.2;
                    norm_arr(i, j, k, dir) = -1.0;
                } else if (idx == 8) {
                    // High
                    vof_arr(i, j, k) = 0.7;
                    tgvof_arr(i, j, k) = 0.6;
                    norm_arr(i, j, k, dir) = 1.0;
                }
            });
    });
}

void init_velocity_etc(amr_wind::Field& velocity)
{
    run_algorithm(velocity, [&](const int lev, const amrex::MFIter& mfi) {
        auto vel_arr = velocity(lev).array(mfi);
        const auto& bx = mfi.validbox();
        amrex::ParallelFor(
            grow(bx, 1), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                // Insert velocity distribution
            });
    });
}

void calc_alpha_flux(
    amr_wind::Field& flux,
    amr_wind::Field& vof,
    amr_wind::Field& tg_vof,
    amr_wind::Field& norm,
    const int& dir,
    const amrex::Real& margin)
{
    run_algorithm(flux, [&](const int lev, const amrex::MFIter& mfi) {
        auto f_arr = flux(lev).array(mfi);
        auto vof_arr = vof(lev).array(mfi);
        auto tgvof_arr = tg_vof(lev).array(mfi);
        auto norm_arr = norm(lev).array(mfi);
        const auto& bx = mfi.validbox();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            f_arr(i, j, k) = amr_wind::overset_ops::alpha_flux(
                i, j, k, dir, margin, vof_arr, tgvof_arr, norm_arr);
        });
    });
}

void calc_velocity_face(
    amr_wind::Field& flux,
    amr_wind::Field& vof,
    amr_wind::Field& velocity,
    const int& dir)
{
    run_algorithm(flux, [&](const int lev, const amrex::MFIter& mfi) {
        auto f_arr = flux(lev).array(mfi);
        auto vof_arr = vof(lev).array(mfi);
        auto vel_arr = velocity(lev).array(mfi);
        const auto& bx = mfi.validbox();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            amrex::Real u_f, v_f, w_f;
            amr_wind::overset_ops::velocity_face(
                i, j, k, dir, vof_arr, vel_arr, u_f, v_f, w_f);
            f_arr(i, j, k, 0) = u_f;
            f_arr(i, j, k, 1) = v_f;
            f_arr(i, j, k, 2) = w_f;
        });
    });
}

amrex::Real check_alpha_flux_impl(amr_wind::Field& flux, const int& dir)
{
    amrex::Real error_total = 0;

    for (int lev = 0; lev < flux.repo().num_active_levels(); ++lev) {

        error_total += amrex::ReduceSum(
            flux(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& f_arr) -> amrex::Real {
                amrex::Real error = 0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                    const int idx = (dir == 0 ? i : (dir == 1 ? j : k));
                    // Difference between actual and expected
                    if (idx == 0) {
                        // Both within margin, will average
                        const amrex::Real flux_answer =
                            0.5 * ((0.5 - 0.41) * -1.0 + (0.5 - 0.55) * -1.0);
                        error += std::abs(f_arr(i, j, k) - flux_answer);
                    } else if (idx == 1) {
                        // Current is low, neighbor within margin
                        const amrex::Real flux_answer = (0.2 - 0.1) * 1.0;
                        error += std::abs(f_arr(i, j, k) - flux_answer);
                    } else if (idx == 2) {
                        // Low on both sides, gphi > 0
                        const amrex::Real flux_answer = (0.2 - 0.1) * -1.0;
                        error += std::abs(f_arr(i, j, k) - flux_answer);
                    } else if (idx == 3) {
                        // Opposite sides of margin, no conditional fits
                        const amrex::Real flux_answer =
                            0.5 * ((0.8 - 0.7) * 1.0 + (0.22 - 0.2) * 1.0);
                        error += std::abs(f_arr(i, j, k) - flux_answer);
                    } else if (idx == 4) {
                        // High on both sides, gphi < 0
                        const amrex::Real flux_answer = (0.8 - 0.7) * -1.0;
                        error += std::abs(f_arr(i, j, k) - flux_answer);
                    } else if (idx == 5) {
                        // High on both sides, gphi > 0
                        const amrex::Real flux_answer = (0.75 - 0.8) * 1.0;
                        error += std::abs(f_arr(i, j, k) - flux_answer);
                    } else if (idx == 6) {
                        // Current is within margin, neighbor is high
                        const amrex::Real flux_answer = (0.75 - 0.8) * -1.0;
                        error += std::abs(f_arr(i, j, k) - flux_answer);
                    } else if (idx == 7) {
                        // Current is low, neighbor is within margin
                        const amrex::Real flux_answer = (0.2 - 0.3) * 1.0;
                        error += std::abs(f_arr(i, j, k) - flux_answer);
                    } else if (idx == 8) {
                        // Opposite sides of margin, no conditional fits
                        const amrex::Real flux_answer =
                            0.5 * ((0.6 - 0.7) * -1.0 + (0.2 - 0.3) * -1.0);
                        error += std::abs(f_arr(i, j, k) - flux_answer);
                    }
                });

                return error;
            });
    }
    return error_total;
}

/*amrex::Real check_velocity_face_impl(amr_wind::Field& flux)
{
    amrex::Real error_total = 0;

    for (int lev = 0; lev < flux.repo().num_active_levels(); ++lev) {

        error_total += amrex::ReduceSum(
            flux(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& f_arr) -> amrex::Real {
                amrex::Real error = 0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                    // Difference between actual and expected
                    const amrex::Real flux_answer = ;
                    error += std::abs(f_arr(i, j, k, 0) - flux_answer);
                    const amrex::Real flux_answer = ;
                    error += std::abs(f_arr(i, j, k, 1) - flux_answer);
                    const amrex::Real flux_answer = ;
                    error += std::abs(f_arr(i, j, k, 2) - flux_answer);
                });

                return error;
            });
    }
    return error_total;
}*/
} // namespace

TEST_F(VOFOversetOps, alpha_flux)
{
    populate_parameters();
    initialize_mesh();

    auto& repo = sim().repo();
    const int ncomp = 1;
    const int nghost = 3;
    auto& vof = repo.declare_field("vof", ncomp, nghost);
    auto& tg_vof = repo.declare_field("target_vof", ncomp, nghost);
    auto& norm = repo.declare_field("int_normal", 3, nghost);
    const amrex::Real margin = 0.1;

    // Create flux fields
    auto& flux_x =
        repo.declare_field("flux_x", 1, 0, 1, amr_wind::FieldLoc::XFACE);
    auto& flux_y =
        repo.declare_field("flux_y", 1, 0, 1, amr_wind::FieldLoc::YFACE);
    auto& flux_z =
        repo.declare_field("flux_z", 1, 0, 1, amr_wind::FieldLoc::ZFACE);

    // -- Variations in x direction -- //
    int dir = 0;
    // Initialize vof and other fields
    init_vof_etc(vof, tg_vof, norm, dir);
    // Populate flux field
    calc_alpha_flux(flux_x, vof, tg_vof, norm, dir, margin);
    // Check results
    amrex::Real error_total = check_alpha_flux_impl(flux_x, dir);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, 1e-15);

    // -- Variations in y direction -- //
    dir = 1;
    // Initialize vof and other fields
    init_vof_etc(vof, tg_vof, norm, dir);
    // Populate flux field
    calc_alpha_flux(flux_y, vof, tg_vof, norm, dir, margin);
    // Check results
    error_total = check_alpha_flux_impl(flux_y, dir);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, 1e-15);

    // -- Variations in z direction -- //
    dir = 2;
    // Initialize vof and other fields
    init_vof_etc(vof, tg_vof, norm, dir);
    // Populate flux field
    calc_alpha_flux(flux_z, vof, tg_vof, norm, dir, margin);
    // Check results
    error_total = check_alpha_flux_impl(flux_z, dir);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, 1e-15);
}

TEST_F(VOFOversetOps, velocity_face) {}

} // namespace amr_wind_tests