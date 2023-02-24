#include "amr-wind/incflo.H"

using namespace amrex;

//
// Print maximum values (useful for tracking evolution)
//
void incflo::PrintMaxValues(const std::string& header)
{
    BL_PROFILE("amr-wind::incflo::PrintMaxValues");

    amrex::Print() << "\nL-inf norm summary: " << header << std::endl
                   << "........................................................"
                      "......................";

    for (int lev = 0; lev <= finest_level; lev++) {
        amrex::Print() << "\nLevel " << lev << std::endl;

        const auto& vel = icns().fields().field;
        amrex::Print() << "  " << std::setw(16) << std::left << vel.name();
        for (int i = 0; i < vel.num_comp(); ++i) {
            amrex::Print() << std::setw(20) << std::right << vel(lev).norm0(i);
        }
        amrex::Print() << std::endl;

        const auto& gradp = grad_p();
        amrex::Print() << "  " << std::setw(16) << std::left << gradp.name();
        for (int i = 0; i < gradp.num_comp(); ++i) {
            amrex::Print() << std::setw(20) << std::right
                           << gradp(lev).norm0(i);
        }
        amrex::Print() << std::endl;

        for (auto& eqn : scalar_eqns()) {
            auto& field = eqn->fields().field;
            amrex::Print() << "  " << std::setw(16) << std::left
                           << field.name();
            for (int i = 0; i < field.num_comp(); ++i) {
                amrex::Print()
                    << std::setw(20) << std::right << field(lev).norm0(i);
            }
            amrex::Print() << std::endl;
        }
    }
    amrex::Print() << "........................................................"
                      "......................"
                   << std::endl
                   << std::endl;
}

void incflo::PrintMaxVelLocations(const std::string& header)
{
    BL_PROFILE("amr-wind::incflo::PrintMaxVelLocations");

    // Get fields
    auto& vel = repo().get_field("velocity");

    // Get infinity norm of velocities
    amrex::Real u_max{-1e8}, v_max{-1e8}, w_max{-1e8};
    // Minima will be negated later
    amrex::Real u_min{-1e8}, v_min{-1e8}, w_min{-1e8};
    for (int lev = 0; lev <= finest_level; lev++) {
        // Use level_mask to only count finest level present
        amrex::iMultiFab level_mask;
        if (lev < finest_level) {
            level_mask = makeFineMask(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                m_sim.mesh().boxArray(lev + 1), amrex::IntVect(2), 1, 0);
        } else {
            level_mask.define(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                1, 0, amrex::MFInfo());
            level_mask.setVal(1);
        }

        u_max = amrex::max(
            u_max,
            amrex::ReduceMax(
                vel(lev), level_mask, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& vel_arr,
                    amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                    amrex::Real max_fab = -1e8;
                    amrex::Loop(
                        bx, [=, &max_fab](int i, int j, int k) noexcept {
                            max_fab = amrex::max(max_fab, vel_arr(i, j, k, 0)) *
                                      (mask_arr(i, j, k) > 0 ? 1.0 : -1.0);
                        });
                    return max_fab;
                }));

        u_min = amrex::max(
            u_min,
            amrex::ReduceMax(
                vel(lev), level_mask, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& vel_arr,
                    amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                    amrex::Real max_fab = -1e8;
                    amrex::Loop(
                        bx, [=, &max_fab](int i, int j, int k) noexcept {
                            max_fab =
                                amrex::max(max_fab, -vel_arr(i, j, k, 0)) *
                                (mask_arr(i, j, k) > 0 ? 1.0 : -1.0);
                        });
                    return max_fab;
                }));

        v_max = amrex::max(
            v_max,
            amrex::ReduceMax(
                vel(lev), level_mask, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& vel_arr,
                    amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                    amrex::Real max_fab = -1e8;
                    amrex::Loop(
                        bx, [=, &max_fab](int i, int j, int k) noexcept {
                            max_fab = amrex::max(max_fab, vel_arr(i, j, k, 1)) *
                                      (mask_arr(i, j, k) > 0 ? 1.0 : -1.0);
                        });
                    return max_fab;
                }));

        v_min = amrex::max(
            v_min,
            amrex::ReduceMax(
                vel(lev), level_mask, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& vel_arr,
                    amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                    amrex::Real max_fab = -1e8;
                    amrex::Loop(
                        bx, [=, &max_fab](int i, int j, int k) noexcept {
                            max_fab =
                                amrex::max(max_fab, -vel_arr(i, j, k, 1)) *
                                (mask_arr(i, j, k) > 0 ? 1.0 : -1.0);
                        });
                    return max_fab;
                }));

        w_max = amrex::max(
            w_max,
            amrex::ReduceMax(
                vel(lev), level_mask, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& vel_arr,
                    amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                    amrex::Real max_fab = -1e8;
                    amrex::Loop(
                        bx, [=, &max_fab](int i, int j, int k) noexcept {
                            max_fab = amrex::max(max_fab, vel_arr(i, j, k, 2)) *
                                      (mask_arr(i, j, k) > 0 ? 1.0 : -1.0);
                        });
                    return max_fab;
                }));

        w_min = amrex::max(
            w_min,
            amrex::ReduceMax(
                vel(lev), level_mask, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& vel_arr,
                    amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                    amrex::Real max_fab = -1e8;
                    amrex::Loop(
                        bx, [=, &max_fab](int i, int j, int k) noexcept {
                            max_fab =
                                amrex::max(max_fab, -vel_arr(i, j, k, 2)) *
                                (mask_arr(i, j, k) > 0 ? 1.0 : -1.0);
                        });
                    return max_fab;
                }));
    }

    // Do additional parallelism stuff
    amrex::ParallelDescriptor::ReduceRealMax(u_max);
    amrex::ParallelDescriptor::ReduceRealMax(v_max);
    amrex::ParallelDescriptor::ReduceRealMax(w_max);
    amrex::ParallelDescriptor::ReduceRealMax(u_min);
    amrex::ParallelDescriptor::ReduceRealMax(v_min);
    amrex::ParallelDescriptor::ReduceRealMax(w_min);

    // Negate minima
    u_min *= -1.0;
    v_min *= -1.0;
    w_min *= -1.0;

    // Get locations of these extrema
    auto problo = (m_sim.mesh().Geom())[0].ProbLoArray();
    auto dx = (m_sim.mesh().Geom())[0].CellSizeArray();
    amrex::GpuArray<amrex::Real, 3> u_max_loc{problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> v_max_loc{problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> w_max_loc{problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> u_min_loc{problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> v_min_loc{problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> w_min_loc{problo[0], problo[1], problo[2]};
    for (int lev = 0; lev <= finest_level; lev++) {
        // Use level_mask to only count finest level present
        amrex::iMultiFab level_mask;
        if (lev < finest_level) {
            level_mask = makeFineMask(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                m_sim.mesh().boxArray(lev + 1), amrex::IntVect(2), 1, 0);
        } else {
            level_mask.define(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                1, 0, amrex::MFInfo());
            level_mask.setVal(1);
        }

        // Loop coordinates directions
        for (int n = 0; n < 3; n++) {
            u_max_loc[n] = amrex::max(
                u_max_loc[n],
                amrex::ReduceMax(
                    vel(lev), level_mask, 0,
                    [=] AMREX_GPU_HOST_DEVICE(
                        amrex::Box const& bx,
                        amrex::Array4<amrex::Real const> const& vel_arr,
                        amrex::Array4<int const> const& mask_arr)
                        -> amrex::Real {
                        amrex::Real loc_fab = problo[n];
                        amrex::Loop(
                            bx, [=, &loc_fab](int i, int j, int k) noexcept {
                                int idx = (n == 0 ? i : (n == 1 ? j : k));
                                amrex::Real offset = 0.5;
                                amrex::Real loc =
                                    problo[n] + (idx + offset) * dx[n];
                                bool mask_check = (mask_arr(i, j, k) > 0);
                                bool loc_check =
                                    (amrex::Math::abs(
                                         u_max - vel_arr(i, j, k, 0)) < 1e-10);
                                loc_fab = amrex::max(
                                    loc_fab,
                                    (mask_check && loc_check ? loc
                                                             : problo[n]));
                            });
                        return loc_fab;
                    }));

            u_min_loc[n] = amrex::max(
                u_min_loc[n],
                amrex::ReduceMax(
                    vel(lev), level_mask, 0,
                    [=] AMREX_GPU_HOST_DEVICE(
                        amrex::Box const& bx,
                        amrex::Array4<amrex::Real const> const& vel_arr,
                        amrex::Array4<int const> const& mask_arr)
                        -> amrex::Real {
                        amrex::Real loc_fab = problo[n];
                        amrex::Loop(
                            bx, [=, &loc_fab](int i, int j, int k) noexcept {
                                int idx = (n == 0 ? i : (n == 1 ? j : k));
                                amrex::Real offset = 0.5;
                                amrex::Real loc =
                                    problo[n] + (idx + offset) * dx[n];
                                bool mask_check = (mask_arr(i, j, k) > 0);
                                bool loc_check =
                                    (amrex::Math::abs(
                                         u_min - vel_arr(i, j, k, 0)) < 1e-10);
                                loc_fab = amrex::max(
                                    loc_fab,
                                    (mask_check && loc_check ? loc
                                                             : problo[n]));
                            });
                        return loc_fab;
                    }));

            v_max_loc[n] = amrex::max(
                v_max_loc[n],
                amrex::ReduceMax(
                    vel(lev), level_mask, 0,
                    [=] AMREX_GPU_HOST_DEVICE(
                        amrex::Box const& bx,
                        amrex::Array4<amrex::Real const> const& vel_arr,
                        amrex::Array4<int const> const& mask_arr)
                        -> amrex::Real {
                        amrex::Real loc_fab = problo[n];
                        amrex::Loop(
                            bx, [=, &loc_fab](int i, int j, int k) noexcept {
                                int idx = (n == 0 ? i : (n == 1 ? j : k));
                                amrex::Real offset = 0.5;
                                amrex::Real loc =
                                    problo[n] + (idx + offset) * dx[n];
                                bool mask_check = (mask_arr(i, j, k) > 0);
                                bool loc_check =
                                    (amrex::Math::abs(
                                         v_max - vel_arr(i, j, k, 1)) < 1e-10);
                                loc_fab = amrex::max(
                                    loc_fab,
                                    (mask_check && loc_check ? loc
                                                             : problo[n]));
                            });
                        return loc_fab;
                    }));

            v_min_loc[n] = amrex::max(
                v_min_loc[n],
                amrex::ReduceMax(
                    vel(lev), level_mask, 0,
                    [=] AMREX_GPU_HOST_DEVICE(
                        amrex::Box const& bx,
                        amrex::Array4<amrex::Real const> const& vel_arr,
                        amrex::Array4<int const> const& mask_arr)
                        -> amrex::Real {
                        amrex::Real loc_fab = problo[n];
                        amrex::Loop(
                            bx, [=, &loc_fab](int i, int j, int k) noexcept {
                                int idx = (n == 0 ? i : (n == 1 ? j : k));
                                amrex::Real offset = 0.5;
                                amrex::Real loc =
                                    problo[n] + (idx + offset) * dx[n];
                                bool mask_check = (mask_arr(i, j, k) > 0);
                                bool loc_check =
                                    (amrex::Math::abs(
                                         v_min - vel_arr(i, j, k, 1)) < 1e-10);
                                loc_fab = amrex::max(
                                    loc_fab,
                                    (mask_check && loc_check ? loc
                                                             : problo[n]));
                            });
                        return loc_fab;
                    }));

            w_max_loc[n] = amrex::max(
                w_max_loc[n],
                amrex::ReduceMax(
                    vel(lev), level_mask, 0,
                    [=] AMREX_GPU_HOST_DEVICE(
                        amrex::Box const& bx,
                        amrex::Array4<amrex::Real const> const& vel_arr,
                        amrex::Array4<int const> const& mask_arr)
                        -> amrex::Real {
                        amrex::Real loc_fab = problo[n];
                        amrex::Loop(
                            bx, [=, &loc_fab](int i, int j, int k) noexcept {
                                int idx = (n == 0 ? i : (n == 1 ? j : k));
                                amrex::Real offset = 0.5;
                                amrex::Real loc =
                                    problo[n] + (idx + offset) * dx[n];
                                bool mask_check = (mask_arr(i, j, k) > 0);
                                bool loc_check =
                                    (amrex::Math::abs(
                                         w_max - vel_arr(i, j, k, 2)) < 1e-10);
                                loc_fab = amrex::max(
                                    loc_fab,
                                    (mask_check && loc_check ? loc
                                                             : problo[n]));
                            });
                        return loc_fab;
                    }));

            w_min_loc[n] = amrex::max(
                w_min_loc[n],
                amrex::ReduceMax(
                    vel(lev), level_mask, 0,
                    [=] AMREX_GPU_HOST_DEVICE(
                        amrex::Box const& bx,
                        amrex::Array4<amrex::Real const> const& vel_arr,
                        amrex::Array4<int const> const& mask_arr)
                        -> amrex::Real {
                        amrex::Real loc_fab = problo[n];
                        amrex::Loop(
                            bx, [=, &loc_fab](int i, int j, int k) noexcept {
                                int idx = (n == 0 ? i : (n == 1 ? j : k));
                                amrex::Real offset = 0.5;
                                amrex::Real loc =
                                    problo[n] + (idx + offset) * dx[n];
                                bool mask_check = (mask_arr(i, j, k) > 0);
                                bool loc_check =
                                    (amrex::Math::abs(
                                         w_min - vel_arr(i, j, k, 2)) < 1e-10);
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
        amrex::ParallelDescriptor::ReduceRealMax(u_max_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(v_max_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(w_max_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(u_min_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(v_min_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(w_min_loc[n]);
    }

    // Output results
    amrex::Print() << "\nL-inf norm vels: " << header << std::endl
                   << "........................................................"
                      "......................" << std::endl;

    amrex::Print() << "Max u: " << std::setw(20) << std::right << u_max;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << u_max_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << u_max_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << u_max_loc[2] << std::endl;
    amrex::Print() << "Min u: " << std::setw(20) << std::right << u_min;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << u_min_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << u_min_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << u_min_loc[2] << std::endl;

    amrex::Print() << "Max v: " << std::setw(20) << std::right << v_max;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << v_max_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << v_max_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << v_max_loc[2] << std::endl;
    amrex::Print() << "Min v: " << std::setw(20) << std::right << v_min;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << v_min_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << v_min_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << v_min_loc[2] << std::endl;

    amrex::Print() << "Max w: " << std::setw(20) << std::right << w_max;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << w_max_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << w_max_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << w_max_loc[2] << std::endl;
    amrex::Print() << "Min w: " << std::setw(20) << std::right << w_min;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << w_min_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << w_min_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << w_min_loc[2] << std::endl;

    amrex::Print() << "........................................................"
                      "......................"
                   << std::endl
                   << std::endl;
}

//
// Print the maximum values of the velocity components and velocity divergence
//
void incflo::PrintMaxVel(int lev) const
{
    BL_PROFILE("amr-wind::incflo::PrintMaxVel");
    amrex::Print() << "max(abs(u/v/w))  = " << velocity()(lev).norm0(0) << "  "
                   << velocity()(lev).norm0(1) << "  "
                   << velocity()(lev).norm0(2) << "  " << std::endl;
}

//
// Print the maximum values of the pressure gradient components and pressure
//
void incflo::PrintMaxGp(int lev) const
{
    BL_PROFILE("amr-wind::incflo::PrintMaxGp");
    amrex::Print() << "max(abs(gpx/gpy/gpz/p))  = " << grad_p()(lev).norm0(0)
                   << "  " << grad_p()(lev).norm0(1) << "  "
                   << grad_p()(lev).norm0(2) << "  " << pressure()(lev).norm0(0)
                   << "  " << std::endl;
}

void incflo::CheckForNans(int lev) const
{
    BL_PROFILE("amr-wind::incflo::CheckForNans");
    bool ro_has_nans = density()(lev).contains_nan(false);
    bool ug_has_nans = velocity()(lev).contains_nan(false);
    bool vg_has_nans = velocity()(lev).contains_nan(true);
    bool wg_has_nans = velocity()(lev).contains_nan(true);
    bool pg_has_nans = pressure()(lev).contains_nan(false);

    if (ro_has_nans) {
        amrex::Print() << "WARNING: ro contains NaNs!!!";
    }

    if (ug_has_nans) {
        amrex::Print() << "WARNING: u contains NaNs!!!";
    }

    if (vg_has_nans) {
        amrex::Print() << "WARNING: v contains NaNs!!!";
    }

    if (wg_has_nans) {
        amrex::Print() << "WARNING: w contains NaNs!!!";
    }

    if (pg_has_nans) {
        amrex::Print() << "WARNING: p contains NaNs!!!";
    }
}
