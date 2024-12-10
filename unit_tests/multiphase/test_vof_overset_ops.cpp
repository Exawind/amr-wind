#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/overset/overset_ops_routines.H"

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

void init_velocity_etc(
    amr_wind::Field& velocity, amr_wind::Field& vof, const int& dir)
{
    run_algorithm(velocity, [&](const int lev, const amrex::MFIter& mfi) {
        auto vel_arr = velocity(lev).array(mfi);
        auto vof_arr = vof(lev).array(mfi);
        const auto& bx = mfi.validbox();
        amrex::ParallelFor(
            grow(bx, 2), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                vel_arr(i, j, k) = 0.0;
                vof_arr(i, j, k) = 0.0;
                const int idx = (dir == 0 ? i : (dir == 1 ? j : k));
                if (idx == -1) {
                    vof_arr(i, j, k) = 0.41;
                    vel_arr(i, j, k, 0) = 1.0;
                    vel_arr(i, j, k, 1) = 2.0;
                    vel_arr(i, j, k, 2) = 3.0;
                } else if (idx == 0) {
                    vof_arr(i, j, k) = 0.55;
                    vel_arr(i, j, k, 0) = 1.5;
                    vel_arr(i, j, k, 1) = 2.5;
                    vel_arr(i, j, k, 2) = 3.5;
                } else if (idx == 1) {
                    vof_arr(i, j, k) = 0.2;
                    vel_arr(i, j, k, 0) = 2.0;
                    vel_arr(i, j, k, 1) = 3.0;
                    vel_arr(i, j, k, 2) = 4.0;
                } else if (idx == 2) {
                    vof_arr(i, j, k) = 0.2;
                    vel_arr(i, j, k, 0) = 2.5;
                    vel_arr(i, j, k, 1) = 3.5;
                    vel_arr(i, j, k, 2) = 4.5;
                }
            });
    });
}

void init_gp_rho_etc(
    amr_wind::Field& gp,
    amr_wind::Field& rho,
    amr_wind::Field& vof,
    const int& dir,
    const amrex::Real& rho1,
    const amrex::Real& rho2)
{
    run_algorithm(gp, [&](const int lev, const amrex::MFIter& mfi) {
        auto gp_arr = gp(lev).array(mfi);
        auto rho_arr = rho(lev).array(mfi);
        auto vof_arr = vof(lev).array(mfi);
        const auto& bx = mfi.validbox();
        amrex::ParallelFor(
            grow(bx, 2), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                gp_arr(i, j, k) = 0.0;
                vof_arr(i, j, k) = 0.0;
                const int idx = (dir == 0 ? i : (dir == 1 ? j : k));
                if (idx == -1) {
                    vof_arr(i, j, k) = 0.41;
                    gp_arr(i, j, k, 0) = 1.0;
                    gp_arr(i, j, k, 1) = 2.0;
                    gp_arr(i, j, k, 2) = 3.0;
                } else if (idx == 0) {
                    vof_arr(i, j, k) = 0.55;
                    gp_arr(i, j, k, 0) = 1.5;
                    gp_arr(i, j, k, 1) = 2.5;
                    gp_arr(i, j, k, 2) = 3.5;
                } else if (idx == 1) {
                    vof_arr(i, j, k) = 0.2;
                    gp_arr(i, j, k, 0) = 2.0;
                    gp_arr(i, j, k, 1) = 3.0;
                    gp_arr(i, j, k, 2) = 4.0;
                } else if (idx == 2) {
                    vof_arr(i, j, k) = 0.2;
                    gp_arr(i, j, k, 0) = 2.5;
                    gp_arr(i, j, k, 1) = 3.5;
                    gp_arr(i, j, k, 2) = 4.5;
                }
                rho_arr(i, j, k) =
                    rho1 * vof_arr(i, j, k) + rho2 * (1. - vof_arr(i, j, k));
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
    amr_wind::Field& velocity,
    amr_wind::Field& vof,
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

void calc_gp_rho_face(
    amr_wind::Field& flux,
    amr_wind::Field& gp,
    amr_wind::Field& rho,
    amr_wind::Field& vof,
    const int& dir,
    const int& ioff = 0,
    bool random_fflag = false)
{
    run_algorithm(flux, [&](const int lev, const amrex::MFIter& mfi) {
        auto f_arr = flux(lev).array(mfi);
        auto vof_arr = vof(lev).array(mfi);
        auto gp_arr = gp(lev).array(mfi);
        auto rho_arr = rho(lev).array(mfi);
        const auto& bx = mfi.validbox();
        amrex::ParallelFor(
            grow(bx, flux.num_grow()),
            [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                amrex::Real u_f, v_f, w_f;
                amrex::Real fac = 1.0;
                if (random_fflag) {
                    fac = ((i + j + k) % 8 == 0) ? 0. : fac;
                    f_arr(i, j, k, 3 + ioff) = fac;
                }
                amr_wind::overset_ops::gp_rho_face(
                    i, j, k, dir, vof_arr, gp_arr, rho_arr, u_f, v_f, w_f);
                f_arr(i, j, k, 0 + ioff) = u_f * fac;
                f_arr(i, j, k, 1 + ioff) = v_f * fac;
                f_arr(i, j, k, 2 + ioff) = w_f * fac;
            });
    });
}

void calc_psrc(
    amr_wind::Field& f_fx,
    amr_wind::Field& f_fy,
    amr_wind::Field& f_fz,
    amr_wind::Field& f_vof,
    amr_wind::Field& f_psrc,
    amr_wind::Field& f_psrc_manual)
{
    constexpr amrex::Real tiny = std::numeric_limits<amrex::Real>::epsilon();
    run_algorithm(f_vof, [&](const int lev, const amrex::MFIter& mfi) {
        auto fx = f_fx(lev).const_array(mfi);
        auto fy = f_fy(lev).const_array(mfi);
        auto fz = f_fz(lev).const_array(mfi);
        auto vof = f_vof(lev).const_array(mfi);
        auto psrc = f_psrc(lev).array(mfi);
        auto psmn = f_psrc_manual(lev).array(mfi);
        const auto& bx = mfi.validbox();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            psrc(i, j, k) = amr_wind::overset_ops::gp_flux_tensor(
                                i, j, k, fx, fy, fz, tiny) &&
                            amr_wind::overset_ops::normal_reinit_tensor(
                                i, j, k, fx, fy, fz, vof, tiny);

            const amrex::Real fxgpx_flux =
                (fx(i, j, k, 5) + fx(i, j, k - 1, 5) + fx(i, j - 1, k, 5) +
                 fx(i, j - 1, k - 1, 5)) /
                (fx(i, j, k, 8) + fx(i, j, k - 1, 8) + fx(i, j - 1, k, 8) +
                 fx(i, j - 1, k - 1, 8) + tiny);
            const amrex::Real fxgpy_flux =
                (fx(i, j, k, 6) + fx(i, j, k - 1, 6) + fx(i, j - 1, k, 6) +
                 fx(i, j - 1, k - 1, 6)) /
                (fx(i, j, k, 8) + fx(i, j, k - 1, 8) + fx(i, j - 1, k, 8) +
                 fx(i, j - 1, k - 1, 8) + tiny);
            const amrex::Real fxgpz_flux =
                (fx(i, j, k, 7) + fx(i, j, k - 1, 7) + fx(i, j - 1, k, 7) +
                 fx(i, j - 1, k - 1, 7)) /
                (fx(i, j, k, 8) + fx(i, j, k - 1, 8) + fx(i, j - 1, k, 8) +
                 fx(i, j - 1, k - 1, 8) + tiny);

            const amrex::Real fygpx_flux =
                (fy(i, j, k, 5) + fy(i - 1, j, k, 5) + fy(i, j, k - 1, 5) +
                 fy(i - 1, j, k - 1, 5)) /
                (fy(i, j, k, 8) + fy(i - 1, j, k, 8) + fy(i, j, k - 1, 8) +
                 fy(i - 1, j, k - 1, 8) + tiny);
            const amrex::Real fygpy_flux =
                (fy(i, j, k, 6) + fy(i - 1, j, k, 6) + fy(i, j, k - 1, 6) +
                 fy(i - 1, j, k - 1, 6)) /
                (fy(i, j, k, 8) + fy(i - 1, j, k, 8) + fy(i, j, k - 1, 8) +
                 fy(i - 1, j, k - 1, 8) + tiny);
            const amrex::Real fygpz_flux =
                (fy(i, j, k, 7) + fy(i - 1, j, k, 7) + fy(i, j, k - 1, 7) +
                 fy(i - 1, j, k - 1, 7)) /
                (fy(i, j, k, 8) + fy(i - 1, j, k, 8) + fy(i, j, k - 1, 8) +
                 fy(i - 1, j, k - 1, 8) + tiny);

            const amrex::Real fzgpx_flux =
                (fz(i, j, k, 5) + fz(i - 1, j, k, 5) + fz(i, j - 1, k, 5) +
                 fz(i - 1, j - 1, k, 5)) /
                (fz(i, j, k, 8) + fz(i - 1, j, k, 8) + fz(i, j - 1, k, 8) +
                 fz(i - 1, j - 1, k, 8) + tiny);
            const amrex::Real fzgpy_flux =
                (fz(i, j, k, 6) + fz(i - 1, j, k, 6) + fz(i, j - 1, k, 6) +
                 fz(i - 1, j - 1, k, 6)) /
                (fz(i, j, k, 8) + fz(i - 1, j, k, 8) + fz(i, j - 1, k, 8) +
                 fz(i - 1, j - 1, k, 8) + tiny);
            const amrex::Real fzgpz_flux =
                (fz(i, j, k, 7) + fz(i - 1, j, k, 7) + fz(i, j - 1, k, 7) +
                 fz(i - 1, j - 1, k, 7)) /
                (fz(i, j, k, 8) + fz(i - 1, j, k, 8) + fz(i, j - 1, k, 8) +
                 fz(i - 1, j - 1, k, 8) + tiny);

            // Calculate interface normal at pressure node
            amrex::Real n0 =
                ((vof(i, j, k) - vof(i - 1, j, k)) * fx(i, j, k, 8) +
                 (vof(i, j - 1, k) - vof(i - 1, j - 1, k)) *
                     fx(i, j - 1, k, 8) +
                 (vof(i, j, k - 1) - vof(i - 1, j, k - 1)) *
                     fx(i, j, k - 1, 8) +
                 (vof(i, j - 1, k - 1) - vof(i - 1, j - 1, k - 1)) *
                     fx(i, j - 1, k - 1, 8)) /
                (fx(i, j, k, 8) + fx(i, j, k - 1, 8) + fx(i, j - 1, k, 8) +
                 fx(i, j - 1, k - 1, 8) + tiny);
            amrex::Real n1 =
                ((vof(i, j, k) - vof(i, j - 1, k)) * fy(i, j, k, 8) +
                 (vof(i - 1, j, k) - vof(i - 1, j - 1, k)) *
                     fy(i - 1, j, k, 8) +
                 (vof(i, j, k - 1) - vof(i, j - 1, k - 1)) *
                     fy(i, j, k - 1, 8) +
                 (vof(i - 1, j, k - 1) - vof(i - 1, j - 1, k - 1)) *
                     fy(i - 1, j, k - 1, 8)) /
                (fy(i, j, k, 8) + fy(i, j, k - 1, 8) + fy(i - 1, j, k, 8) +
                 fy(i - 1, j, k - 1, 8) + tiny);
            amrex::Real n2 =
                ((vof(i, j, k) - vof(i, j, k - 1)) * fz(i, j, k, 8) +
                 (vof(i - 1, j, k) - vof(i - 1, j, k - 1)) *
                     fz(i - 1, j, k, 8) +
                 (vof(i, j - 1, k) - vof(i, j - 1, k - 1)) *
                     fz(i, j - 1, k, 8) +
                 (vof(i - 1, j - 1, k) - vof(i - 1, j - 1, k - 1)) *
                     fz(i - 1, j - 1, k, 8)) /
                (fz(i, j, k, 8) + fz(i, j - 1, k, 8) + fz(i - 1, j, k, 8) +
                 fz(i - 1, j - 1, k, 8) + tiny);
            // Normalize vector
            amrex::Real nmag = std::sqrt(n0 * n0 + n1 * n1 + n2 * n2) + tiny;
            n0 /= nmag;
            n1 /= nmag;
            n2 /= nmag;

            // Perform double contraction
            psmn(i, j, k) = fxgpx_flux * n0 * n0 + fxgpy_flux * n0 * n1 +
                            fxgpz_flux * n0 * n2 + fygpx_flux * n1 * n0 +
                            fygpy_flux * n1 * n1 + fygpz_flux * n1 * n2 +
                            fzgpx_flux * n2 * n0 + fzgpy_flux * n2 * n1 +
                            fzgpz_flux * n2 * n2;
        });
    });
}

amrex::Real check_alpha_flux_impl(amr_wind::Field& flux, const int& dir)
{
    amrex::Real error_total = 0;

    for (int lev = 0; lev < flux.repo().num_active_levels(); ++lev) {

        error_total += amrex::ReduceSum(
            flux(lev), 0,
            [=] AMREX_GPU_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& f_arr) -> amrex::Real {
                amrex::Real error = 0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                    const int idx = (dir == 0 ? i : (dir == 1 ? j : k));
                    // Difference between actual and expected
                    if (idx == 0) {
                        // Both within margin, will average
                        const amrex::Real flux_answer =
                            0.5 * ((0.5 - 0.41) * 1.0 + (0.5 - 0.55) * -1.0);
                        error += std::abs(f_arr(i, j, k) - flux_answer);
                    } else if (idx == 1) {
                        // Current is low, neighbor within margin
                        const amrex::Real flux_answer = (0.2 - 0.1) * -1.0;
                        error += std::abs(f_arr(i, j, k) - flux_answer);
                    } else if (idx == 2) {
                        // Low on both sides, gphi > 0
                        const amrex::Real flux_answer = (0.2 - 0.1) * 1.0;
                        error += std::abs(f_arr(i, j, k) - flux_answer);
                    } else if (idx == 3) {
                        // Opposite sides of margin, no conditional fits
                        const amrex::Real flux_answer =
                            0.5 * ((0.8 - 0.7) * -1.0 + (0.22 - 0.2) * 1.0);
                        error += std::abs(f_arr(i, j, k) - flux_answer);
                    } else if (idx == 4) {
                        // High on both sides, gphi < 0
                        const amrex::Real flux_answer = (0.8 - 0.7) * 1.0;
                        error += std::abs(f_arr(i, j, k) - flux_answer);
                    } else if (idx == 5) {
                        // High on both sides, gphi > 0
                        const amrex::Real flux_answer = (0.75 - 0.8) * -1.0;
                        error += std::abs(f_arr(i, j, k) - flux_answer);
                    } else if (idx == 6) {
                        // Current is within margin, neighbor is high
                        const amrex::Real flux_answer = (0.75 - 0.8) * 1.0;
                        error += std::abs(f_arr(i, j, k) - flux_answer);
                    } else if (idx == 7) {
                        // Current is low, neighbor is within margin
                        const amrex::Real flux_answer = (0.2 - 0.3) * -1.0;
                        error += std::abs(f_arr(i, j, k) - flux_answer);
                    } else if (idx == 8) {
                        // Opposite sides of margin, no conditional fits
                        const amrex::Real flux_answer =
                            0.5 * ((0.6 - 0.7) * -1.0 + (0.2 - 0.3) * 1.0);
                        error += std::abs(f_arr(i, j, k) - flux_answer);
                    }
                });

                return error;
            });
    }
    return error_total;
}

amrex::Real check_velocity_face_impl(amr_wind::Field& flux, const int& dir)
{
    amrex::Real error_total = 0;

    for (int lev = 0; lev < flux.repo().num_active_levels(); ++lev) {

        error_total += amrex::ReduceSum(
            flux(lev), 0,
            [=] AMREX_GPU_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& f_arr) -> amrex::Real {
                amrex::Real error = 0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                    const int idx = (dir == 0 ? i : (dir == 1 ? j : k));
                    if (idx == 0) {
                        // gphi > 0, uwpind from the "left"
                        amrex::Real flux_answer = 1.0;
                        error += std::abs(f_arr(i, j, k, 0) - flux_answer);
                        flux_answer = 2.0;
                        error += std::abs(f_arr(i, j, k, 1) - flux_answer);
                        flux_answer = 3.0;
                        error += std::abs(f_arr(i, j, k, 2) - flux_answer);
                    } else if (idx == 1) {
                        // gphi < 0, upwind from the "right"
                        amrex::Real flux_answer = 2.0;
                        error += std::abs(f_arr(i, j, k, 0) - flux_answer);
                        flux_answer = 3.0;
                        error += std::abs(f_arr(i, j, k, 1) - flux_answer);
                        flux_answer = 4.0;
                        error += std::abs(f_arr(i, j, k, 2) - flux_answer);
                    } else if (idx == 2) {
                        // gphi = 0, average both sides
                        amrex::Real flux_answer = 0.5 * (2.5 + 2.0);
                        error += std::abs(f_arr(i, j, k, 0) - flux_answer);
                        flux_answer = 0.5 * (3.5 + 3.0);
                        error += std::abs(f_arr(i, j, k, 1) - flux_answer);
                        flux_answer = 0.5 * (4.5 + 4.0);
                        error += std::abs(f_arr(i, j, k, 2) - flux_answer);
                    }
                });

                return error;
            });
    }
    return error_total;
}

amrex::Real check_gp_rho_face_impl(
    amr_wind::Field& flux,
    const int& dir,
    const amrex::Real& rho1,
    const amrex::Real& rho2)
{
    amrex::Real error_total = 0;

    for (int lev = 0; lev < flux.repo().num_active_levels(); ++lev) {

        error_total += amrex::ReduceSum(
            flux(lev), 0,
            [=] AMREX_GPU_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& f_arr) -> amrex::Real {
                amrex::Real error = 0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                    const int idx = (dir == 0 ? i : (dir == 1 ? j : k));
                    if (idx == 0) {
                        // gphi > 0, uwpind from the "left"
                        const amrex::Real rho_answer =
                            rho1 * 0.41 + rho2 * (1. - 0.41);
                        amrex::Real flux_answer = 1.0 / rho_answer;
                        error += std::abs(f_arr(i, j, k, 0) - flux_answer);
                        flux_answer = 2.0 / rho_answer;
                        error += std::abs(f_arr(i, j, k, 1) - flux_answer);
                        flux_answer = 3.0 / rho_answer;
                        error += std::abs(f_arr(i, j, k, 2) - flux_answer);
                    } else if (idx == 1) {
                        // gphi < 0, upwind from the "right"
                        const amrex::Real rho_answer =
                            rho1 * 0.2 + rho2 * (1. - 0.2);
                        amrex::Real flux_answer = 2.0 / rho_answer;
                        error += std::abs(f_arr(i, j, k, 0) - flux_answer);
                        flux_answer = 3.0 / rho_answer;
                        error += std::abs(f_arr(i, j, k, 1) - flux_answer);
                        flux_answer = 4.0 / rho_answer;
                        error += std::abs(f_arr(i, j, k, 2) - flux_answer);
                    } else if (idx == 2) {
                        // gphi = 0, average both sides
                        const amrex::Real rho_r =
                            rho1 * 0.2 + rho2 * (1. - 0.2);
                        const amrex::Real rho_l =
                            rho1 * 0.2 + rho2 * (1. - 0.2);
                        amrex::Real flux_answer =
                            0.5 * (2.5 / rho_r + 2.0 / rho_l);
                        error += std::abs(f_arr(i, j, k, 0) - flux_answer);
                        flux_answer = 0.5 * (3.5 / rho_r + 3.0 / rho_l);
                        error += std::abs(f_arr(i, j, k, 1) - flux_answer);
                        flux_answer = 0.5 * (4.5 / rho_r + 4.0 / rho_l);
                        error += std::abs(f_arr(i, j, k, 2) - flux_answer);
                    }
                });

                return error;
            });
    }
    return error_total;
}

amrex::Real
check_psrc_manual_impl(amr_wind::Field& psrc, amr_wind::Field& psrc_manual)
{
    amrex::Real error_total = 0;

    for (int lev = 0; lev < psrc.repo().num_active_levels(); ++lev) {

        error_total += amrex::ReduceSum(
            psrc(lev), psrc_manual(lev), 0,
            [=] AMREX_GPU_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& f_arr,
                amrex::Array4<amrex::Real const> const& fm_arr) -> amrex::Real {
                amrex::Real error = 0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                    error += std::abs(f_arr(i, j, k) - fm_arr(i, j, k));
                });

                return error;
            });
    }
    return error_total;
}
} // namespace

TEST_F(VOFOversetOps, alpha_flux)
{
    populate_parameters();
    initialize_mesh();

    auto& repo = sim().repo();
    const int nghost = 3;
    auto& vof = repo.declare_field("vof", 1, nghost);
    auto& tg_vof = repo.declare_field("target_vof", 1, nghost);
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

TEST_F(VOFOversetOps, velocity_face)
{
    populate_parameters();
    initialize_mesh();

    auto& repo = sim().repo();
    const int nghost = 3;
    auto& velocity = repo.declare_field("velocity", 3, nghost);
    auto& vof = repo.declare_field("vof", 1, nghost);

    // Create flux fields
    auto& flux_x =
        repo.declare_field("flux_x", 3, 0, 1, amr_wind::FieldLoc::XFACE);
    auto& flux_y =
        repo.declare_field("flux_y", 3, 0, 1, amr_wind::FieldLoc::YFACE);
    auto& flux_z =
        repo.declare_field("flux_z", 3, 0, 1, amr_wind::FieldLoc::ZFACE);

    // -- Variations in x direction -- //
    int dir = 0;
    // Initialize velocity and vof
    init_velocity_etc(velocity, vof, dir);
    // Populate flux field
    calc_velocity_face(flux_x, velocity, vof, dir);
    // Check results
    amrex::Real error_total = check_velocity_face_impl(flux_x, dir);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, 1e-15);

    // -- Variations in y direction -- //
    dir = 1;
    // Initialize velocity and vof
    init_velocity_etc(velocity, vof, dir);
    // Populate flux field
    calc_velocity_face(flux_y, velocity, vof, dir);
    // Check results
    error_total = check_velocity_face_impl(flux_y, dir);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, 1e-15);

    // -- Variations in z direction -- //
    dir = 2;
    // Initialize velocity and vof
    init_velocity_etc(velocity, vof, dir);
    // Populate flux field
    calc_velocity_face(flux_z, velocity, vof, dir);
    // Check results
    error_total = check_velocity_face_impl(flux_z, dir);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, 1e-15);
}

TEST_F(VOFOversetOps, gp_rho_face)
{
    populate_parameters();
    initialize_mesh();

    const amrex::Real rho_liq = 1000.;
    const amrex::Real rho_gas = 1.;

    auto& repo = sim().repo();
    const int nghost = 3;
    auto& gp = repo.declare_field("gp", 3, nghost);
    auto& rho = repo.declare_field("density", 1, nghost);
    auto& vof = repo.declare_field("vof", 1, nghost);

    // Create flux fields
    auto& flux_x =
        repo.declare_field("flux_x", 3, 0, 1, amr_wind::FieldLoc::XFACE);
    auto& flux_y =
        repo.declare_field("flux_y", 3, 0, 1, amr_wind::FieldLoc::YFACE);
    auto& flux_z =
        repo.declare_field("flux_z", 3, 0, 1, amr_wind::FieldLoc::ZFACE);

    // -- Variations in x direction -- //
    int dir = 0;
    // Initialize gp, rho, and vof
    init_gp_rho_etc(gp, rho, vof, dir, rho_liq, rho_gas);
    // Populate flux field
    calc_gp_rho_face(flux_x, gp, rho, vof, dir);
    // Check results
    amrex::Real error_total =
        check_gp_rho_face_impl(flux_x, dir, rho_liq, rho_gas);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, 1e-15);

    // -- Variations in y direction -- //
    dir = 1;
    // Initialize gp, rho, and vof
    init_gp_rho_etc(gp, rho, vof, dir, rho_liq, rho_gas);
    // Populate flux field
    calc_gp_rho_face(flux_y, gp, rho, vof, dir);
    // Check results
    error_total = check_gp_rho_face_impl(flux_y, dir, rho_liq, rho_gas);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, 1e-15);

    // -- Variations in z direction -- //
    dir = 2;
    // Initialize gp, rho, and vof
    init_gp_rho_etc(gp, rho, vof, dir, rho_liq, rho_gas);
    // Populate flux field
    calc_gp_rho_face(flux_z, gp, rho, vof, dir);
    // Check results
    error_total = check_gp_rho_face_impl(flux_z, dir, rho_liq, rho_gas);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, 1e-15);
}

TEST_F(VOFOversetOps, pseudo_vscale_dt)
{
    populate_parameters();
    initialize_mesh();

    auto& repo = sim().repo();
    const int nghost = 3;
    auto& vof = repo.declare_field("vof", 1, nghost);
    auto& tg_vof = repo.declare_field("target_vof", 1, nghost);
    auto& norm = repo.declare_field("int_normal", 3, nghost);
    auto& iblank = repo.declare_int_field("iblank_cell", 1, nghost);
    iblank.setVal(-1);
    constexpr amrex::Real margin = 0.1;
    constexpr amrex::Real convg_tol = 1e-8;
    // With vof and target_vof arrays, max vof removed is 50%, doubling pdt
    constexpr amrex::Real pdt_answer = 2.0;
    // With a single level, pseudo velocity scale should be dx of lev 0
    const auto dx_lev0 = repo.mesh().Geom(0).CellSizeArray();
    const amrex::Real pvs_answer =
        std::min(std::min(dx_lev0[0], dx_lev0[1]), dx_lev0[2]);

    // Pseudo-velocity scale, should be the smallest dx in iblank region
    const auto& iblank_cell = repo.get_int_field("iblank_cell");
    const amrex::Real max_pvscale = 100.;
    amrex::Real pvscale = max_pvscale;
    for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
        const amrex::Real pvscale_lev =
            amr_wind::overset_ops::calculate_pseudo_velocity_scale(
                iblank_cell(lev), repo.mesh().Geom(lev).CellSizeArray(),
                max_pvscale);
        pvscale = std::min(pvscale, pvscale_lev);
    }
    amrex::ParallelDescriptor::ReduceRealMin(pvscale);
    EXPECT_DOUBLE_EQ(pvscale, pvs_answer);

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
    // Calculate pseudo dt
    amrex::Real ptfac = 100.0;
    for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
        const amrex::Real ptfac_lev =
            amr_wind::overset_ops::calculate_pseudo_dt_flux(
                flux_x(lev), flux_y(lev), flux_z(lev), vof(lev), dx_lev0,
                convg_tol) /
            pvscale;
        ptfac = amrex::min(ptfac, ptfac_lev);
    }
    amrex::ParallelDescriptor::ReduceRealMin(ptfac);
    EXPECT_DOUBLE_EQ(ptfac, pdt_answer);
    // Zero flux for subsequent tests
    flux_x.setVal(0.0);

    // -- Variations in y direction -- //
    dir = 1;
    // Initialize vof and other fields
    init_vof_etc(vof, tg_vof, norm, dir);
    // Populate flux field
    calc_alpha_flux(flux_y, vof, tg_vof, norm, dir, margin);
    // Calculate pseudo dt
    ptfac = 100.0;
    for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
        const amrex::Real ptfac_lev =
            amr_wind::overset_ops::calculate_pseudo_dt_flux(
                flux_x(lev), flux_y(lev), flux_z(lev), vof(lev), dx_lev0,
                convg_tol) /
            pvscale;
        ptfac = amrex::min(ptfac, ptfac_lev);
    }
    amrex::ParallelDescriptor::ReduceRealMin(ptfac);
    EXPECT_DOUBLE_EQ(ptfac, pdt_answer);
    // Zero flux for subsequent test
    flux_y.setVal(0.0);

    // -- Variations in z direction -- //
    dir = 2;
    // Initialize vof and other fields
    init_vof_etc(vof, tg_vof, norm, dir);
    // Populate flux field
    calc_alpha_flux(flux_z, vof, tg_vof, norm, dir, margin);
    // Calculate pseudo dt
    ptfac = 100.0;
    for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
        const amrex::Real ptfac_lev =
            amr_wind::overset_ops::calculate_pseudo_dt_flux(
                flux_x(lev), flux_y(lev), flux_z(lev), vof(lev), dx_lev0,
                convg_tol) /
            pvscale;
        ptfac = amrex::min(ptfac, ptfac_lev);
    }
    amrex::ParallelDescriptor::ReduceRealMin(ptfac);
    EXPECT_DOUBLE_EQ(ptfac, pdt_answer);

    // Redo z direction with max being 1, as in overall algorithm
    ptfac = 1.0;
    for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
        const amrex::Real ptfac_lev =
            amr_wind::overset_ops::calculate_pseudo_dt_flux(
                flux_x(lev), flux_y(lev), flux_z(lev), vof(lev), dx_lev0,
                convg_tol) /
            pvscale;
        ptfac = amrex::min(ptfac, ptfac_lev);
    }
    amrex::ParallelDescriptor::ReduceRealMin(ptfac);
    EXPECT_DOUBLE_EQ(ptfac, 1.0);
}

TEST_F(VOFOversetOps, psource_manual)
{
    populate_parameters();
    initialize_mesh();

    const amrex::Real rho_liq = 1000.;
    const amrex::Real rho_gas = 1.;

    auto& repo = sim().repo();
    const int nghost = 3;
    auto& gp = repo.declare_field("gp", 3, nghost);
    auto& rho = repo.declare_field("density", 1, nghost);
    auto& vof = repo.declare_field("vof", 1, nghost);

    auto& psrc = repo.declare_field("psrc", 1, 0, 1, amr_wind::FieldLoc::NODE);
    auto& psmn =
        repo.declare_field("psrc_manual", 1, 0, 1, amr_wind::FieldLoc::NODE);

    // Create flux fields
    auto& flux_x =
        repo.declare_field("flux_x", 9, 1, 1, amr_wind::FieldLoc::XFACE);
    auto& flux_y =
        repo.declare_field("flux_y", 9, 1, 1, amr_wind::FieldLoc::YFACE);
    auto& flux_z =
        repo.declare_field("flux_z", 9, 1, 1, amr_wind::FieldLoc::ZFACE);

    // Initialize gp_rho fluxes in every direction
    int dir = 0;
    init_gp_rho_etc(gp, rho, vof, dir, rho_liq, rho_gas);
    calc_gp_rho_face(flux_x, gp, rho, vof, dir, 5, true);
    dir = 1;
    init_gp_rho_etc(gp, rho, vof, dir, rho_liq, rho_gas);
    calc_gp_rho_face(flux_y, gp, rho, vof, dir, 5, true);
    dir = 2;
    init_gp_rho_etc(gp, rho, vof, dir, rho_liq, rho_gas);
    calc_gp_rho_face(flux_z, gp, rho, vof, dir, 5, true);

    // Calculate psrc with vector/tensor ops and manually
    calc_psrc(flux_x, flux_y, flux_z, vof, psrc, psmn);

    // Check difference
    amrex::Real error_total = check_psrc_manual_impl(psrc, psmn);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, 1e-10);
}

} // namespace amr_wind_tests