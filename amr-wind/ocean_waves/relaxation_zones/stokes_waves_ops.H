#ifndef STOKES_WAVES_OPS_H
#define STOKES_WAVES_OPS_H

#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/ocean_waves/relaxation_zones/StokesWaves.H"
#include "amr-wind/ocean_waves/relaxation_zones/stokes_waves_K.H"
#include "amr-wind/ocean_waves/OceanWavesOps.H"
#include "amr-wind/ocean_waves/OceanWaves.H"
#include "amr-wind/ocean_waves/relaxation_zones/relaxation_zones_ops.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"

namespace amr_wind::ocean_waves::ops {

template <>
struct ReadInputsOp<StokesWaves>
{
    void operator()(
        StokesWaves::DataType& data, const ::amr_wind::utils::MultiParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();
        relaxation_zones::read_inputs(wdata, info, pp);

        // Get gravity from MultiPhase physics, assume negative z
        wdata.g = -data.sim().physics_manager().get<MultiPhase>().gravity()[2];
        // Get wave attributes
        pp.get("wave_height", wdata.wave_height);
        pp.get("order", wdata.order);
        if (pp.contains("wave_length")) {
            pp.get("wave_length", wdata.wave_length);
        } else {
            // Wave length will be calculated using wave period
            amrex::Real wave_period = 1.0;
            pp.get("wave_period", wave_period);
            // User can specify a different order for wavelength calculation,
            // e.g., an order of 1 will use the linear dispersion relation.
            // Default is to use the same order as the waves themselves.
            int order_lambda = wdata.order;
            pp.query("stokes_wavelength_order", order_lambda);
            // Get user-specified convergence criteria, if supplied
            amrex::Real tol_lambda = 1e-10;
            int itmax_lambda = 40;
            pp.query("stokes_wavelength_tolerance", tol_lambda);
            pp.query("stokes_wavelength_iter_max", itmax_lambda);
            wdata.wave_length = relaxation_zones::stokes_wave_length(
                wave_period, wdata.water_depth, wdata.wave_height, order_lambda,
                wdata.g, tol_lambda, itmax_lambda);
            // Abort if wave length is negative
            if (wdata.wave_length < 0.0) {
                if (wdata.wave_length == -1) {
                    amrex::Print() << "Stokes wave length too close to 0.\n";
                }
                if (wdata.wave_length == -2) {
                    amrex::Print() << "Stokes wave length is not a number.\n";
                }
                if (wdata.wave_length == -3) {
                    amrex::Print() << "Stokes wave length calculation used "
                                      "maximum iterations before converging.\n";
                }
                amrex::Abort(
                    "Failed to properly calculate a wave length in "
                    "stokes_waves_ops.H");
            } else {
                amrex::Print() << "Stokes wave length calculated to be "
                               << wdata.wave_length << " m.\n";
            }
        }
    }
};

template <>
struct InitDataOp<StokesWaves>
{
    void operator()(
        StokesWaves::DataType& data, int level, const amrex::Geometry& geom)
    {
        const auto& wdata = data.meta();

        auto& sim = data.sim();
        // cppcheck-suppress constVariableReference
        auto& m_levelset = sim.repo().get_field("levelset");
        // cppcheck-suppress constVariableReference
        auto& velocity = sim.repo().get_field("velocity");
        const auto& problo = geom.ProbLoArray();
        const auto& probhi = geom.ProbHiArray();
        const auto& dx = geom.CellSizeArray();
        for (amrex::MFIter mfi(m_levelset(level)); mfi.isValid(); ++mfi) {

            const auto& phi = m_levelset(level).array(mfi);
            const auto& vel = velocity(level).array(mfi);

            const amrex::Real wave_height = wdata.wave_height;
            const amrex::Real wave_length = wdata.wave_length;
            const amrex::Real water_depth = wdata.water_depth;
            const amrex::Real zero_sea_level = wdata.zsl;
            const amrex::Real gen_length = wdata.gen_length;
            const amrex::Real beach_length = wdata.beach_length;
            const amrex::Real g = wdata.g;
            const int order = wdata.order;
            const bool has_beach = wdata.has_beach;
            const bool init_wave_field = wdata.init_wave_field;
            const auto& gbx3 = mfi.growntilebox(3);

            amrex::ParallelFor(
                gbx3, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                    // Wave profile
                    amrex::Real eta_w{0.0}, u_w{0.0}, v_w{0.0}, w_w{0.0};
                    relaxation_zones::stokes_waves(
                        order, wave_length, water_depth, wave_height,
                        zero_sea_level, g, x, z, 0.0, eta_w, u_w, v_w, w_w);
                    const utils::WaveVec wave_sol{u_w, v_w, w_w, eta_w};

                    // Quiescent profile
                    const utils::WaveVec quiescent{
                        0.0, 0.0, 0.0, zero_sea_level};

                    // Specify initial state for each region of domain
                    const auto bulk = init_wave_field ? wave_sol : quiescent;
                    const auto outlet = has_beach ? quiescent : wave_sol;

                    auto local_profile = utils::harmonize_profiles_1d(
                        x, problo[0], gen_length, probhi[0], beach_length,
                        wave_sol, bulk, outlet);

                    phi(i, j, k) = local_profile[3] - z;
                    const amrex::Real cell_length_2D =
                        std::sqrt(dx[0] * dx[0] + dx[2] * dx[2]);
                    if (phi(i, j, k) + cell_length_2D >= 0) {
                        vel(i, j, k, 0) = local_profile[0];
                        vel(i, j, k, 1) = local_profile[1];
                        vel(i, j, k, 2) = local_profile[2];
                    }
                });
        }
    }
};

template <>
struct UpdateRelaxZonesOp<StokesWaves>
{
    void operator()(StokesWaves::DataType& data)
    {
        const auto& wdata = data.meta();

        auto& sim = data.sim();

        const auto& time = sim.time().new_time();
        // cppcheck-suppress constVariableReference
        auto& m_ow_levelset = sim.repo().get_field("ow_levelset");
        // cppcheck-suppress constVariableReference
        auto& m_ow_velocity = sim.repo().get_field("ow_velocity");

        auto nlevels = sim.repo().num_active_levels();
        auto geom = sim.mesh().Geom();

        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& problo = geom[lev].ProbLoArray();
            const auto& dx = geom[lev].CellSizeArray();

            for (amrex::MFIter mfi(m_ow_levelset(lev)); mfi.isValid(); ++mfi) {
                const auto& phi = m_ow_levelset(lev).array(mfi);
                const auto& vel = m_ow_velocity(lev).array(mfi);

                const amrex::Real wave_height = wdata.wave_height;
                const amrex::Real wave_length = wdata.wave_length;
                const amrex::Real water_depth = wdata.water_depth;
                const amrex::Real zero_sea_level = wdata.zsl;
                const amrex::Real g = wdata.g;
                const int order = wdata.order;

                const auto& gbx = mfi.growntilebox();
                amrex::ParallelFor(
                    gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                        amrex::Real eta{0.0}, u_w{0.0}, v_w{0.0}, w_w{0.0};

                        relaxation_zones::stokes_waves(
                            order, wave_length, water_depth, wave_height,
                            zero_sea_level, g, x, z, time, eta, u_w, v_w, w_w);

                        phi(i, j, k) = eta - z;
                        const amrex::Real cell_length_2D =
                            std::sqrt(dx[0] * dx[0] + dx[2] * dx[2]);
                        if (phi(i, j, k) + cell_length_2D >= 0) {
                            // Wave velocity within a cell of interface
                            vel(i, j, k, 0) = u_w;
                            vel(i, j, k, 1) = v_w;
                            vel(i, j, k, 2) = w_w;
                        } else {
                            vel(i, j, k, 0) = 0.;
                            vel(i, j, k, 1) = 0.;
                            vel(i, j, k, 2) = 0.;
                        }
                    });
            }
        }
    }
};

} // namespace amr_wind::ocean_waves::ops

#endif /* STOKES_WAVES_OPS_H */
