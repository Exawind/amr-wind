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

        pp.get("wave_length", wdata.wave_length);
        pp.get("wave_height", wdata.wave_height);
        pp.get("order", wdata.order);
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
        // cppcheck-suppress constVariable
        auto& m_levelset = sim.repo().get_field("levelset");
        // cppcheck-suppress constVariable
        auto& velocity = sim.repo().get_field("velocity");
        const auto& problo = geom.ProbLoArray();
        const auto& dx = geom.CellSizeArray();
        for (amrex::MFIter mfi(m_levelset(level)); mfi.isValid(); ++mfi) {

            auto phi = m_levelset(level).array(mfi);
            auto vel = velocity(level).array(mfi);

            const amrex::Real zsl = wdata.zsl;
            const auto& gbx3 = mfi.growntilebox(3);
            amrex::ParallelFor(
                gbx3, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                    phi(i, j, k) = zsl - z;
                    vel(i, j, k, 0) = 0.0;
                    vel(i, j, k, 1) = 0.0;
                    vel(i, j, k, 2) = 0.0;
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
        // cppcheck-suppress constVariable
        auto& m_ow_levelset = sim.repo().get_field("ow_levelset");
        // cppcheck-suppress constVariable
        auto& m_ow_velocity = sim.repo().get_field("ow_velocity");

        auto nlevels = sim.repo().num_active_levels();
        auto geom = sim.mesh().Geom();

        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& problo = geom[lev].ProbLoArray();
            const auto& dx = geom[lev].CellSizeArray();

            for (amrex::MFIter mfi(m_ow_levelset(lev)); mfi.isValid(); ++mfi) {
                auto phi = m_ow_levelset(lev).array(mfi);
                auto vel = m_ow_velocity(lev).array(mfi);

                const amrex::Real waveheight = wdata.wave_height;
                const amrex::Real wavelength = wdata.wave_length;
                const amrex::Real waterdepth = wdata.water_depth;
                const int order = wdata.order;

                const auto& gbx = mfi.growntilebox();
                amrex::ParallelFor(
                    gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                        amrex::Real eta{0.0}, u_w{0.0}, v_w{0.0}, w_w{0.0};

                        relaxation_zones::stokes_waves(
                            order, wavelength, waterdepth, waveheight, x, z,
                            time, eta, u_w, v_w, w_w);

                        phi(i, j, k) = eta - z;
                        if (phi(i, j, k) + 0.5 * dx[2] >= 0) {
                            vel(i, j, k, 0) = u_w;
                            vel(i, j, k, 1) = v_w;
                            vel(i, j, k, 2) = w_w;
                        }
                    });
            }
        }
    }
};

} // namespace amr_wind::ocean_waves::ops

#endif /* STOKES_WAVES_OPS_H */
