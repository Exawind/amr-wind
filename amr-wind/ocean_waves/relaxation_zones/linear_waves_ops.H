#ifndef LINEAR_WAVES_OPS_H
#define LINEAR_WAVES_OPS_H

#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/ocean_waves/relaxation_zones/LinearWaves.H"
#include "amr-wind/ocean_waves/OceanWavesOps.H"
#include "amr-wind/ocean_waves/OceanWaves.H"
#include "amr-wind/ocean_waves/relaxation_zones/relaxation_zones_ops.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"

#include <cmath>

#include "AMReX_ParmParse.H"
#include "AMReX_Print.H"
#include "AMReX_ParallelDescriptor.H"

namespace amr_wind::ocean_waves::ops {

template <>
struct ReadInputsOp<LinearWaves>
{
    void operator()(
        LinearWaves::DataType& data, const ::amr_wind::utils::MultiParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();
        relaxation_zones::read_inputs(wdata, info, pp);

        pp.get("wave_length", wdata.wave_length);
        pp.get("wave_height", wdata.wave_height);
    }
};

template <>
struct InitDataOp<LinearWaves>
{
    void operator()(
        LinearWaves::DataType& data, int level, const amrex::Geometry& geom)
    {
        const auto& wdata = data.meta();

        auto& sim = data.sim();

        auto& m_levelset = sim.repo().get_field("levelset");
        // cppcheck-suppress constVariable
        auto& m_velocity = sim.repo().get_field("velocity");

        const auto& problo = geom.ProbLoArray();
        const auto& dx = geom.CellSizeArray();
        for (amrex::MFIter mfi(m_levelset(level)); mfi.isValid(); ++mfi) {

            auto phi = m_levelset(level).array(mfi);
            auto vel = m_velocity(level).array(mfi);

            const amrex::Real zsl = wdata.zsl;
            const auto& gbx3 = mfi.growntilebox(3);

            if (wdata.init_wave_field) {
                const amrex::Real waveheight = wdata.wave_height;
                const amrex::Real wavelength = wdata.wave_length;
                const amrex::Real waterdepth = wdata.water_depth;
                amrex::ParallelFor(
                    gbx3, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real xc = problo[0] + (i + 0.5) * dx[0];
                        const amrex::Real zc = problo[2] + (k + 0.5) * dx[2];

                        const amrex::Real wavenumber = 2. * M_PI / wavelength;
                        const amrex::Real omega = std::pow(
                            wavenumber * 9.81 *
                                std::tanh(wavenumber * waterdepth),
                            0.5);
                        const amrex::Real phase = wavenumber * xc;
                        const amrex::Real eta =
                            waveheight / 2.0 * std::cos(phase);

                        phi(i, j, k) = eta - zc;

                        if (phi(i, j, k) >= 0) {
                            vel(i, j, k, 0) =
                                omega * waveheight / 2.0 *
                                std::cosh(wavenumber * (zc + waterdepth)) /
                                std::sinh(wavenumber * waterdepth) *
                                std::cos(phase);
                            vel(i, j, k, 1) = 0.0;
                            vel(i, j, k, 2) =
                                omega * waveheight / 2.0 *
                                std::sinh(wavenumber * (zc + waterdepth)) /
                                std::sinh(wavenumber * waterdepth) *
                                std::sin(phase);
                        }
                    });

            } else {
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
    }
};

template <>
struct UpdateRelaxZonesOp<LinearWaves>
{
    void operator()(LinearWaves::DataType& data)
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
                const amrex::Real grav = wdata.grav;

                const auto& gbx = mfi.growntilebox(3);
                amrex::ParallelFor(
                    gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real xc = problo[0] + (i + 0.5) * dx[0];
                        const amrex::Real zc = problo[2] + (k + 0.5) * dx[2];

                        const amrex::Real wavenumber = 2. * M_PI / wavelength;
                        const amrex::Real omega = std::pow(
                            wavenumber * grav *
                                std::tanh(wavenumber * waterdepth),
                            0.5);
                        const amrex::Real phase =
                            wavenumber * xc - omega * time;

                        const amrex::Real eta =
                            waveheight / 2.0 * std::cos(phase);

                        phi(i, j, k) = eta - zc;

                        if (phi(i, j, k) + 0.5 * dx[2] >= 0) {
                            vel(i, j, k, 0) =
                                omega * waveheight / 2.0 *
                                std::cosh(wavenumber * (zc + waterdepth)) /
                                std::sinh(wavenumber * waterdepth) *
                                std::cos(phase);
                            vel(i, j, k, 1) = 0.0;
                            vel(i, j, k, 2) =
                                omega * waveheight / 2.0 *
                                std::sinh(wavenumber * (zc + waterdepth)) /
                                std::sinh(wavenumber * waterdepth) *
                                std::sin(phase);
                        }
                    });
            }
        }
    }
};

template <>
struct ApplyPressureBCOp<LinearWaves>
{
    void operator()(LinearWaves::DataType& data)
    {

        constexpr int idim = 0;

        const auto& wdata = data.meta();
        auto& sim = data.sim();
        // Getting the time at n+1/2
        const auto& time = sim.time().new_time();
        auto& pressure = sim.repo().get_field("p");

        const int nlevels = sim.repo().num_active_levels();

        amrex::Orientation xlo(amrex::Direction::x, amrex::Orientation::low);
        amrex::Orientation xhi(amrex::Direction::x, amrex::Orientation::high);
        if ((pressure.bc_type()[xlo] != BC::wave_generation) &&
            (pressure.bc_type()[xhi] != BC::wave_generation)) {
            return;
        }

        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& geom = sim.mesh().Geom(lev);
            const auto& domain = geom.Domain();
            const auto& problo = geom.ProbLoArray();
            const auto& dx = geom.CellSizeArray();
            amrex::MFItInfo mfi_info{};

            auto& pres_lev = pressure(lev);

            if (amrex::Gpu::notInLaunchRegion()) {
                mfi_info.SetDynamic(true);
            }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(pres_lev, mfi_info); mfi.isValid(); ++mfi) {
                const amrex::Real waveheight = wdata.wave_height;
                const amrex::Real wavelength = wdata.wave_length;
                const amrex::Real waterdepth = wdata.water_depth;
                const amrex::Real grav = wdata.grav;
                const amrex::Real rho = wdata.water_density;

                const auto& nbx = mfi.nodaltilebox();
                auto parr = pres_lev.array(mfi);
                // if (nbx.smallEnd(idim) == domain.smallEnd(idim) &&
                //     pressure.bc_type()[xlo] == BC::wave_generation) {
                amrex::ParallelFor(
                    nbx,
                    // amrex::bdryLo(nbx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real xn = problo[0] + i * dx[0];
                        const amrex::Real zn = problo[2] + k * dx[2];

                        const amrex::Real wavenumber = 2. * M_PI / wavelength;
                        const amrex::Real omega = std::pow(
                            wavenumber * grav *
                                std::tanh(wavenumber * waterdepth),
                            0.5);
                        const amrex::Real phase =
                            wavenumber * xn - omega * time;
                        const amrex::Real eta =
                            waveheight / 2.0 * std::cos(phase);
                        if (eta - zn + 0.5 * dx[2] >= 0) {
                            parr(i, j, k) =
                                grav *
                                (waterdepth +
                                 waveheight / 2. *
                                     std::cosh(wavenumber * (zn + waterdepth)) /
                                     std::cosh(wavenumber * waterdepth) *
                                     std::cos(phase));
                        } else {
                            parr(i, j, k) = 0.;
                        }
                    });
                //}
                /*
                if (nbx.bigEnd(idim) == domain.bigEnd(idim) &&
                    pressure.bc_type()[xhi] == BC::wave_generation) {
                    amrex::ParallelFor(
                        amrex::bdryHi(nbx, idim),
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            const amrex::Real zn = problo[2] + k * dx[2];

                            if (zn < problo[2] + waterdepth) {
                                parr(i, j, k) =
                                    waterdepth * grav * rho; // Static pressure
                            }
                        });
                */
            }
        }
    }
};

} // namespace amr_wind::ocean_waves::ops

#endif /* LINEAR_WAVES_OPS_H */
