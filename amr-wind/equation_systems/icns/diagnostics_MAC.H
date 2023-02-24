#ifndef DIAGNOSTICS_MAC_H
#define DIAGNOSTICS_MAC_H

#include "AMReX_MultiFabUtil.H"

namespace amr_wind {

static void PrintMaxMACVelLocations(const amr_wind::FieldRepo& repo, const std::string& header)
{
    BL_PROFILE("amr-wind::incflo::PrintMaxMACVelLocations");

    // Get fields
    auto& u_mac = repo.get_field("u_mac");
    auto& v_mac = repo.get_field("v_mac");
    auto& w_mac = repo.get_field("w_mac");
    const int finest_level = repo.num_active_levels() - 1;

    // Get infinity norm of mac velocities
    amrex::Real uMAC_max{-1e8}, vMAC_max{-1e8}, wMAC_max{-1e8};
    // Minima will be negated later
    amrex::Real uMAC_min{-1e8}, vMAC_min{-1e8}, wMAC_min{-1e8};
    for (int lev = 0; lev <= finest_level; lev++) {
        // Use level_mask to only count finest level present
        amrex::iMultiFab level_mask;
        if (lev < finest_level) {
            level_mask = makeFineMask(
                repo.mesh().boxArray(lev), repo.mesh().DistributionMap(lev),
                repo.mesh().boxArray(lev + 1), amrex::IntVect(2), 1, 0);
        } else {
            level_mask.define(
                repo.mesh().boxArray(lev), repo.mesh().DistributionMap(lev),
                1, 0, amrex::MFInfo());
            level_mask.setVal(1);
        }

        uMAC_max = amrex::max(
            uMAC_max,
            amrex::ReduceMax(
                u_mac(lev), level_mask, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& um_arr,
                    amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                    amrex::Real max_fab = -1e8;
                    amrex::Loop(
                        bx, [=, &max_fab](int i, int j, int k) noexcept {
                            max_fab =
                                amrex::max(max_fab, um_arr(i, j, k)) *
                                (mask_arr(i, j, k) + mask_arr(i - 1, j, k) > 0
                                     ? 1.0
                                     : -1.0);
                        });
                    return max_fab;
                }));

        uMAC_min = amrex::max(
            uMAC_min,
            amrex::ReduceMax(
                u_mac(lev), level_mask, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& um_arr,
                    amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                    amrex::Real max_fab = -1e8;
                    amrex::Loop(
                        bx, [=, &max_fab](int i, int j, int k) noexcept {
                            max_fab =
                                amrex::max(max_fab, -um_arr(i, j, k)) *
                                (mask_arr(i, j, k) + mask_arr(i - 1, j, k) > 0
                                     ? 1.0
                                     : -1.0);
                        });
                    return max_fab;
                }));

        vMAC_max = amrex::max(
            vMAC_max,
            amrex::ReduceMax(
                v_mac(lev), level_mask, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& vm_arr,
                    amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                    amrex::Real max_fab = -1e8;
                    amrex::Loop(
                        bx, [=, &max_fab](int i, int j, int k) noexcept {
                            max_fab =
                                amrex::max(max_fab, vm_arr(i, j, k)) *
                                (mask_arr(i, j, k) + mask_arr(i, j - 1, k) > 0
                                     ? 1.0
                                     : -1.0);
                        });
                    return max_fab;
                }));

        vMAC_min = amrex::max(
            vMAC_min,
            amrex::ReduceMax(
                v_mac(lev), level_mask, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& vm_arr,
                    amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                    amrex::Real max_fab = -1e8;
                    amrex::Loop(
                        bx, [=, &max_fab](int i, int j, int k) noexcept {
                            max_fab =
                                amrex::max(max_fab, -vm_arr(i, j, k)) *
                                (mask_arr(i, j, k) + mask_arr(i, j - 1, k) > 0
                                     ? 1.0
                                     : -1.0);
                        });
                    return max_fab;
                }));

        wMAC_max = amrex::max(
            wMAC_max,
            amrex::ReduceMax(
                w_mac(lev), level_mask, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& wm_arr,
                    amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                    amrex::Real max_fab = -1e8;
                    amrex::Loop(
                        bx, [=, &max_fab](int i, int j, int k) noexcept {
                            max_fab =
                                amrex::max(max_fab, wm_arr(i, j, k)) *
                                (mask_arr(i, j, k) + mask_arr(i, j, k - 1) > 0
                                     ? 1.0
                                     : -1.0);
                        });
                    return max_fab;
                }));

        wMAC_min = amrex::max(
            wMAC_min,
            amrex::ReduceMax(
                w_mac(lev), level_mask, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& wm_arr,
                    amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                    amrex::Real max_fab = -1e8;
                    amrex::Loop(
                        bx, [=, &max_fab](int i, int j, int k) noexcept {
                            max_fab =
                                amrex::max(max_fab, -wm_arr(i, j, k)) *
                                (mask_arr(i, j, k) + mask_arr(i, j, k - 1) > 0
                                     ? 1.0
                                     : -1.0);
                        });
                    return max_fab;
                }));
    }

    // Do additional parallelism stuff
    amrex::ParallelDescriptor::ReduceRealMax(uMAC_max);
    amrex::ParallelDescriptor::ReduceRealMax(vMAC_max);
    amrex::ParallelDescriptor::ReduceRealMax(wMAC_max);
    amrex::ParallelDescriptor::ReduceRealMax(uMAC_min);
    amrex::ParallelDescriptor::ReduceRealMax(vMAC_min);
    amrex::ParallelDescriptor::ReduceRealMax(wMAC_min);

    // Negate minima
    uMAC_min *= -1.0;
    vMAC_min *= -1.0;
    wMAC_min *= -1.0;

    // Get locations of these extrema
    auto problo = (repo.mesh().Geom())[0].ProbLoArray();
    auto dx = (repo.mesh().Geom())[0].CellSizeArray();
    amrex::GpuArray<amrex::Real, 3> uMAC_max_loc{
        problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> vMAC_max_loc{
        problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> wMAC_max_loc{
        problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> uMAC_min_loc{
        problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> vMAC_min_loc{
        problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> wMAC_min_loc{
        problo[0], problo[1], problo[2]};
    for (int lev = 0; lev <= finest_level; lev++) {
        // Use level_mask to only count finest level present
        amrex::iMultiFab level_mask;
        if (lev < finest_level) {
            level_mask = makeFineMask(
                repo.mesh().boxArray(lev), repo.mesh().DistributionMap(lev),
                repo.mesh().boxArray(lev + 1), amrex::IntVect(2), 1, 0);
        } else {
            level_mask.define(
                repo.mesh().boxArray(lev), repo.mesh().DistributionMap(lev),
                1, 0, amrex::MFInfo());
            level_mask.setVal(1);
        }

        // Loop coordinates directions
        for (int n = 0; n < 3; n++) {
            uMAC_max_loc[n] = amrex::max(
                uMAC_max_loc[n],
                amrex::ReduceMax(
                    u_mac(lev), level_mask, 0,
                    [=] AMREX_GPU_HOST_DEVICE(
                        amrex::Box const& bx,
                        amrex::Array4<amrex::Real const> const& um_arr,
                        amrex::Array4<int const> const& mask_arr)
                        -> amrex::Real {
                        amrex::Real loc_fab = problo[n];
                        amrex::Loop(
                            bx, [=, &loc_fab](int i, int j, int k) noexcept {
                                int idx = (n == 0 ? i : (n == 1 ? j : k));
                                amrex::Real offset = (n == 0 ? 0.0 : 0.5);
                                amrex::Real loc =
                                    problo[n] + (idx + offset) * dx[n];
                                bool mask_check =
                                    (mask_arr(i, j, k) + mask_arr(i - 1, j, k) >
                                     0);
                                bool loc_check =
                                    (amrex::Math::abs(
                                         uMAC_max - um_arr(i, j, k)) < 1e-10);
                                loc_fab = amrex::max(
                                    loc_fab,
                                    (mask_check && loc_check ? loc
                                                             : problo[n]));
                            });
                        return loc_fab;
                    }));

            uMAC_min_loc[n] = amrex::max(
                uMAC_min_loc[n],
                amrex::ReduceMax(
                    u_mac(lev), level_mask, 0,
                    [=] AMREX_GPU_HOST_DEVICE(
                        amrex::Box const& bx,
                        amrex::Array4<amrex::Real const> const& um_arr,
                        amrex::Array4<int const> const& mask_arr)
                        -> amrex::Real {
                        amrex::Real loc_fab = problo[n];
                        amrex::Loop(
                            bx, [=, &loc_fab](int i, int j, int k) noexcept {
                                int idx = (n == 0 ? i : (n == 1 ? j : k));
                                amrex::Real offset = (n == 0 ? 0.0 : 0.5);
                                amrex::Real loc =
                                    problo[n] + (idx + offset) * dx[n];
                                bool mask_check =
                                    (mask_arr(i, j, k) + mask_arr(i - 1, j, k) >
                                     0);
                                bool loc_check =
                                    (amrex::Math::abs(
                                         uMAC_min - um_arr(i, j, k)) < 1e-10);
                                loc_fab = amrex::max(
                                    loc_fab,
                                    (mask_check && loc_check ? loc
                                                             : problo[n]));
                            });
                        return loc_fab;
                    }));

            vMAC_max_loc[n] = amrex::max(
                vMAC_max_loc[n],
                amrex::ReduceMax(
                    v_mac(lev), level_mask, 0,
                    [=] AMREX_GPU_HOST_DEVICE(
                        amrex::Box const& bx,
                        amrex::Array4<amrex::Real const> const& vm_arr,
                        amrex::Array4<int const> const& mask_arr)
                        -> amrex::Real {
                        amrex::Real loc_fab = problo[n];
                        amrex::Loop(
                            bx, [=, &loc_fab](int i, int j, int k) noexcept {
                                int idx = (n == 0 ? i : (n == 1 ? j : k));
                                amrex::Real offset = (n == 1 ? 0.0 : 0.5);
                                amrex::Real loc =
                                    problo[n] + (idx + offset) * dx[n];
                                bool mask_check =
                                    (mask_arr(i, j, k) + mask_arr(i, j - 1, k) >
                                     0);
                                bool loc_check =
                                    (amrex::Math::abs(
                                         vMAC_max - vm_arr(i, j, k)) < 1e-10);
                                loc_fab = amrex::max(
                                    loc_fab,
                                    (mask_check && loc_check ? loc
                                                             : problo[n]));
                            });
                        return loc_fab;
                    }));

            vMAC_min_loc[n] = amrex::max(
                vMAC_min_loc[n],
                amrex::ReduceMax(
                    v_mac(lev), level_mask, 0,
                    [=] AMREX_GPU_HOST_DEVICE(
                        amrex::Box const& bx,
                        amrex::Array4<amrex::Real const> const& vm_arr,
                        amrex::Array4<int const> const& mask_arr)
                        -> amrex::Real {
                        amrex::Real loc_fab = problo[n];
                        amrex::Loop(
                            bx, [=, &loc_fab](int i, int j, int k) noexcept {
                                int idx = (n == 0 ? i : (n == 1 ? j : k));
                                amrex::Real offset = (n == 1 ? 0.0 : 0.5);
                                amrex::Real loc =
                                    problo[n] + (idx + offset) * dx[n];
                                bool mask_check =
                                    (mask_arr(i, j, k) + mask_arr(i, j - 1, k) >
                                     0);
                                bool loc_check =
                                    (amrex::Math::abs(
                                         vMAC_min - vm_arr(i, j, k)) < 1e-10);
                                loc_fab = amrex::max(
                                    loc_fab,
                                    (mask_check && loc_check ? loc
                                                             : problo[n]));
                            });
                        return loc_fab;
                    }));

            wMAC_max_loc[n] = amrex::max(
                wMAC_max_loc[n],
                amrex::ReduceMax(
                    w_mac(lev), level_mask, 0,
                    [=] AMREX_GPU_HOST_DEVICE(
                        amrex::Box const& bx,
                        amrex::Array4<amrex::Real const> const& wm_arr,
                        amrex::Array4<int const> const& mask_arr)
                        -> amrex::Real {
                        amrex::Real loc_fab = problo[n];
                        amrex::Loop(
                            bx, [=, &loc_fab](int i, int j, int k) noexcept {
                                int idx = (n == 0 ? i : (n == 1 ? j : k));
                                amrex::Real offset = (n == 2 ? 0.0 : 0.5);
                                amrex::Real loc =
                                    problo[n] + (idx + offset) * dx[n];
                                bool mask_check =
                                    (mask_arr(i, j, k) + mask_arr(i, j, k - 1) >
                                     0);
                                bool loc_check =
                                    (amrex::Math::abs(
                                         wMAC_max - wm_arr(i, j, k)) < 1e-10);
                                loc_fab = amrex::max(
                                    loc_fab,
                                    (mask_check && loc_check ? loc
                                                             : problo[n]));
                            });
                        return loc_fab;
                    }));

            wMAC_min_loc[n] = amrex::max(
                wMAC_min_loc[n],
                amrex::ReduceMax(
                    w_mac(lev), level_mask, 0,
                    [=] AMREX_GPU_HOST_DEVICE(
                        amrex::Box const& bx,
                        amrex::Array4<amrex::Real const> const& wm_arr,
                        amrex::Array4<int const> const& mask_arr)
                        -> amrex::Real {
                        amrex::Real loc_fab = problo[n];
                        amrex::Loop(
                            bx, [=, &loc_fab](int i, int j, int k) noexcept {
                                int idx = (n == 0 ? i : (n == 1 ? j : k));
                                amrex::Real offset = (n == 2 ? 0.0 : 0.5);
                                amrex::Real loc =
                                    problo[n] + (idx + offset) * dx[n];
                                bool mask_check =
                                    (mask_arr(i, j, k) + mask_arr(i, j, k - 1) >
                                     0);
                                bool loc_check =
                                    (amrex::Math::abs(
                                         wMAC_min - wm_arr(i, j, k)) < 1e-10);
                                loc_fab = amrex::max(
                                    loc_fab,
                                    (mask_check && loc_check ? loc
                                                             : problo[n]));
                            });
                        return loc_fab;
                    }));
        }
    }

    // Additional parallelism
    for (int n = 0; n < 3; ++n) {
        amrex::ParallelDescriptor::ReduceRealMax(uMAC_max_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(vMAC_max_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(wMAC_max_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(uMAC_min_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(vMAC_min_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(wMAC_min_loc[n]);
    }

    // Output results
    amrex::Print() << "\nL-inf norm MAC vels: " << header << std::endl
                   << "........................................................"
                      "......................" << std::endl;

    amrex::Print() << "Max u: " << std::setw(20) << std::right << uMAC_max;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << uMAC_max_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << uMAC_max_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << uMAC_max_loc[2]
                   << std::endl;
    amrex::Print() << "Min u: " << std::setw(20) << std::right << uMAC_min;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << uMAC_min_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << uMAC_min_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << uMAC_min_loc[2]
                   << std::endl;

    amrex::Print() << "Max v: " << std::setw(20) << std::right << vMAC_max;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << vMAC_max_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << vMAC_max_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << vMAC_max_loc[2]
                   << std::endl;
    amrex::Print() << "Min v: " << std::setw(20) << std::right << vMAC_min;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << vMAC_min_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << vMAC_min_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << vMAC_min_loc[2]
                   << std::endl;

    amrex::Print() << "Max w: " << std::setw(20) << std::right << wMAC_max;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << wMAC_max_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << wMAC_max_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << wMAC_max_loc[2]
                   << std::endl;
    amrex::Print() << "Min w: " << std::setw(20) << std::right << wMAC_min;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << wMAC_min_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << wMAC_min_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << wMAC_min_loc[2]
                   << std::endl;

    amrex::Print() << "........................................................"
                      "......................"
                   << std::endl
                   << std::endl;
}

}

#endif