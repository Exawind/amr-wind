#include "amr-wind/ocean_waves/relaxation_zones/relaxation_zones_ops.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "amr-wind/core/MultiParser.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/io_utils.H"

#include "amr-wind/fvm/gradient.H"
#include "amr-wind/core/field_ops.H"

#include "amr-wind/ocean_waves/utils/wave_utils_K.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::ocean_waves::relaxation_zones {

void read_inputs(
    RelaxZonesBaseData& wdata,
    OceanWavesInfo& /*unused*/,
    const ::amr_wind::utils::MultiParser& pp)
{
    // Free surface zero level
    pp.query("zero_sea_level", wdata.zsl);
    pp.query("water_depth", wdata.water_depth);

    // Wave generation/absorption parameters
    pp.query("relax_zone_gen_length", wdata.gen_length);
    if (pp.contains("relax_zone_out_length")) {
        wdata.has_beach = false;
        wdata.has_outprofile = true;
        wdata.init_wave_field = true;
        pp.query("relax_zone_out_length", wdata.beach_length);
    } else {
        pp.query("numerical_beach_length", wdata.beach_length);
        wdata.has_beach = true;
        pp.query("initialize_wave_field", wdata.init_wave_field);
    }

    pp.query("timeramp", wdata.has_ramp);
    if (wdata.has_ramp) {
        pp.query("timeramp_perior", wdata.ramp_period);
    }
}

void init_data_structures(RelaxZonesBaseData& /*unused*/) {}

void apply_relaxation_zones(CFDSim& sim, const RelaxZonesBaseData& wdata)
{
    const int nlevels = sim.repo().num_active_levels();
    auto& m_ow_levelset = sim.repo().get_field("ow_levelset");
    auto& m_ow_vof = sim.repo().get_field("ow_vof");
    const auto& m_ow_vel = sim.repo().get_field("ow_velocity");
    const auto& geom = sim.mesh().Geom();

    const auto& mphase = sim.physics_manager().get<MultiPhase>();
    const amrex::Real rho1 = mphase.rho1();
    const amrex::Real rho2 = mphase.rho2();

    for (int lev = 0; lev < nlevels; ++lev) {
        auto& ls = m_ow_levelset(lev);
        auto& target_vof = m_ow_vof(lev);
        const auto& dx = geom[lev].CellSizeArray();

        for (amrex::MFIter mfi(ls); mfi.isValid(); ++mfi) {
            const auto& gbx = mfi.growntilebox(2);
            const amrex::Array4<amrex::Real>& phi = ls.array(mfi);
            const amrex::Array4<amrex::Real>& volfrac = target_vof.array(mfi);
            const amrex::Real eps = 2. * std::cbrt(dx[0] * dx[1] * dx[2]);
            amrex::ParallelFor(
                gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    volfrac(i, j, k) =
                        multiphase::levelset_to_vof(i, j, k, eps, phi);
                });
        }
    }

    // Get time
    const auto& time = sim.time().new_time();
    const amrex::Real rampf =
        (wdata.has_ramp) ? utils::ramp(time, wdata.ramp_period) : 1.0;

    auto& vof = sim.repo().get_field("vof");
    auto& velocity = sim.repo().get_field("velocity");
    auto& density = sim.repo().get_field("density");

    for (int lev = 0; lev < nlevels; ++lev) {
        for (amrex::MFIter mfi(vof(lev)); mfi.isValid(); ++mfi) {
            const auto& gbx = mfi.growntilebox(2);
            const auto& dx = geom[lev].CellSizeArray();
            const auto& problo = geom[lev].ProbLoArray();
            const auto& probhi = geom[lev].ProbHiArray();
            auto vel = velocity(lev).array(mfi);
            auto rho = density(lev).array(mfi);
            auto volfrac = vof(lev).array(mfi);
            auto target_volfrac = m_ow_vof(lev).array(mfi);
            auto target_vel = m_ow_vel(lev).array(mfi);

            const amrex::Real gen_length = wdata.gen_length;
            const amrex::Real beach_length = wdata.beach_length;
            const amrex::Real zsl = wdata.zsl;
            const bool has_beach = wdata.has_beach;
            const bool has_outprofile = wdata.has_outprofile;

            amrex::ParallelFor(
                gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = amrex::min(
                        amrex::max(problo[0] + (i + 0.5) * dx[0], problo[0]),
                        probhi[0]);
                    const amrex::Real z = amrex::min(
                        amrex::max(problo[2] + (k + 0.5) * dx[2], problo[2]),
                        probhi[2]);

                    // Generation region
                    if (x <= problo[0] + gen_length) {
                        const amrex::Real Gamma =
                            utils::Gamma_generate(x - problo[0], gen_length);
                        const amrex::Real vf =
                            (1. - Gamma) * target_volfrac(i, j, k) * rampf +
                            Gamma * volfrac(i, j, k);
                        volfrac(i, j, k) = (vf > 1. - 1.e-10) ? 1.0 : vf;
                        // Conserve momentum when density changes
                        amrex::Real rho_ = rho1 * volfrac(i, j, k) +
                                           rho2 * (1.0 - volfrac(i, j, k));
                        vel(i, j, k, 0) *= rho(i, j, k) / rho_;
                        vel(i, j, k, 1) *= rho(i, j, k) / rho_;
                        vel(i, j, k, 2) *= rho(i, j, k) / rho_;
                        // Force liquid velocity, update according to mom.
                        vel(i, j, k, 0) =
                            (rho1 * volfrac(i, j, k) *
                                 (rampf * (1. - Gamma) *
                                      target_vel(i, j, k, 0) +
                                  Gamma * vel(i, j, k, 0)) +
                             rho2 * (1. - volfrac(i, j, k)) * vel(i, j, k, 0)) /
                            rho_;
                        vel(i, j, k, 1) =
                            (rho1 * volfrac(i, j, k) *
                                 (rampf * (1. - Gamma) *
                                      target_vel(i, j, k, 1) +
                                  Gamma * vel(i, j, k, 1)) +
                             rho2 * (1. - volfrac(i, j, k)) * vel(i, j, k, 1)) /
                            rho_;
                        vel(i, j, k, 2) =
                            (rho1 * volfrac(i, j, k) *
                                 (rampf * (1. - Gamma) *
                                      target_vel(i, j, k, 2) +
                                  Gamma * vel(i, j, k, 2)) +
                             rho2 * (1. - volfrac(i, j, k)) * vel(i, j, k, 2)) /
                            rho_;
                    }
                    // Numerical beach (sponge layer)
                    if (x + beach_length >= probhi[0]) {
                        const amrex::Real Gamma = utils::Gamma_absorb(
                            x - (probhi[0] - beach_length), beach_length, 1.0);
                        if (has_beach) {
                            volfrac(i, j, k) =
                                (1.0 - Gamma) *
                                    utils::free_surface_to_vof(zsl, z, dx[2]) +
                                Gamma * volfrac(i, j, k);
                            vel(i, j, k, 0) =
                                Gamma * vel(i, j, k, 0) * volfrac(i, j, k);
                            vel(i, j, k, 1) =
                                Gamma * vel(i, j, k, 1) * volfrac(i, j, k);
                            vel(i, j, k, 2) =
                                Gamma * vel(i, j, k, 2) * volfrac(i, j, k);
                        }
                        if (has_outprofile) {
                            const amrex::Real vf =
                                (1. - Gamma) * target_volfrac(i, j, k) * rampf +
                                Gamma * volfrac(i, j, k);
                            volfrac(i, j, k) = (vf > 1. - 1.e-10) ? 1.0 : vf;
                            // Conserve momentum when density changes
                            amrex::Real rho_ = rho1 * volfrac(i, j, k) +
                                               rho2 * (1.0 - volfrac(i, j, k));
                            vel(i, j, k, 0) *= rho(i, j, k) / rho_;
                            vel(i, j, k, 1) *= rho(i, j, k) / rho_;
                            vel(i, j, k, 2) *= rho(i, j, k) / rho_;
                            // Force liquid velocity, update according to mom.
                            vel(i, j, k, 0) = (rho1 * volfrac(i, j, k) *
                                                   (rampf * (1. - Gamma) *
                                                        target_vel(i, j, k, 0) +
                                                    Gamma * vel(i, j, k, 0)) +
                                               rho2 * (1. - volfrac(i, j, k)) *
                                                   vel(i, j, k, 0)) /
                                              rho_;
                            vel(i, j, k, 1) = (rho1 * volfrac(i, j, k) *
                                                   (rampf * (1. - Gamma) *
                                                        target_vel(i, j, k, 1) +
                                                    Gamma * vel(i, j, k, 1)) +
                                               rho2 * (1. - volfrac(i, j, k)) *
                                                   vel(i, j, k, 1)) /
                                              rho_;
                            vel(i, j, k, 2) = (rho1 * volfrac(i, j, k) *
                                                   (rampf * (1. - Gamma) *
                                                        target_vel(i, j, k, 2) +
                                                    Gamma * vel(i, j, k, 2)) +
                                               rho2 * (1. - volfrac(i, j, k)) *
                                                   vel(i, j, k, 2)) /
                                              rho_;
                        }
                    }

                    // Make sure that density is updated before entering the
                    // solution
                    rho(i, j, k) = rho1 * volfrac(i, j, k) +
                                   rho2 * (1. - volfrac(i, j, k));
                });
        }
    }
    // This helps for having periodic boundaries, but will need to be addressed
    // for the general case
    vof.fillpatch(0.0);
    velocity.fillpatch(0.0);
    density.fillpatch(0.0);
}

void prepare_netcdf_file(
    const std::string& ncfile,
    const RelaxZonesBaseData& meta,
    const OceanWavesInfo& info)
{
    amrex::ignore_unused(ncfile, meta, info);
}

void write_netcdf(
    const std::string& ncfile,
    const RelaxZonesBaseData& meta,
    const OceanWavesInfo& info,
    const amrex::Real time)
{
    amrex::ignore_unused(ncfile, meta, info, time);
}

} // namespace amr_wind::ocean_waves::relaxation_zones
