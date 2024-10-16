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
        pp.query("numerical_beach_length_factor", wdata.beach_length_factor);
        pp.query("initialize_wave_field", wdata.init_wave_field);
    }

    wdata.has_ramp = pp.contains("timeramp_period");
    if (wdata.has_ramp) {
        pp.get("timeramp_period", wdata.ramp_period);
    }

    amrex::ParmParse pp_multiphase("MultiPhase");
    pp_multiphase.add("water_level", wdata.zsl);
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
    constexpr amrex::Real vof_tiny = 1e-12;

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
            const amrex::Real beach_length_factor = wdata.beach_length_factor;
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
                            utils::gamma_generate(x - problo[0], gen_length);
                        // Get bounded new vof, incorporate with increment
                        amrex::Real new_vof =
                            (1. - Gamma) * target_volfrac(i, j, k) +
                            Gamma * volfrac(i, j, k);
                        new_vof = (new_vof > 1. - vof_tiny)
                                      ? 1.0
                                      : (new_vof < vof_tiny ? 0.0 : new_vof);
                        const amrex::Real dvf = new_vof - volfrac(i, j, k);
                        volfrac(i, j, k) += rampf * dvf;
                        // Force liquid velocity only if target vof present
                        const amrex::Real fvel_liq =
                            (target_volfrac(i, j, k) > vof_tiny) ? 1.0 : 0.0;
                        amrex::Real rho_ = rho1 * volfrac(i, j, k) +
                                           rho2 * (1.0 - volfrac(i, j, k));
                        for (int n = 0; n < vel.ncomp; ++n) {
                            // Get updated liquid velocity
                            amrex::Real vel_liq = vel(i, j, k, n);
                            const amrex::Real dvel_liq =
                                ((1. - Gamma) * target_vel(i, j, k, n) +
                                 Gamma * vel_liq) -
                                vel_liq;
                            vel_liq += rampf * fvel_liq * dvel_liq;
                            // If liquid was added, that liquid has target_vel
                            amrex::Real integrated_vel_liq =
                                volfrac(i, j, k) * vel_liq;
                            integrated_vel_liq +=
                                rampf * fvel_liq * amrex::max(0.0, dvf) *
                                (target_vel(i, j, k, n) - vel(i, j, k, n));
                            // Update overall velocity using momentum
                            vel(i, j, k, n) = (rho1 * integrated_vel_liq +
                                               rho2 * (1. - volfrac(i, j, k)) *
                                                   vel(i, j, k, n)) /
                                              rho_;
                        }
                    }
                    // Outlet region
                    if (x + beach_length >= probhi[0]) {
                        const amrex::Real Gamma = utils::gamma_absorb(
                            x - (probhi[0] - beach_length), beach_length,
                            beach_length_factor);
                        // Numerical beach (sponge layer)
                        if (has_beach) {
                            // Get bounded new vof, save increment
                            amrex::Real new_vof =
                                (1. - Gamma) *
                                    utils::free_surface_to_vof(zsl, z, dx[2]) +
                                Gamma * volfrac(i, j, k);
                            new_vof =
                                (new_vof > 1. - vof_tiny)
                                    ? 1.0
                                    : (new_vof < vof_tiny ? 0.0 : new_vof);
                            const amrex::Real dvf = new_vof - volfrac(i, j, k);
                            volfrac(i, j, k) = new_vof;
                            // Conserve momentum when density changes
                            amrex::Real rho_ = rho1 * volfrac(i, j, k) +
                                               rho2 * (1.0 - volfrac(i, j, k));
                            // Target solution in liquid is vel = 0, assume
                            // added liquid already has target velocity
                            for (int n = 0; n < vel.ncomp; ++n) {
                                vel(i, j, k, n) =
                                    (rho1 * (volfrac(i, j, k) * Gamma -
                                             amrex::max(0.0, dvf)) +
                                     rho2 * (1. - volfrac(i, j, k))) *
                                    vel(i, j, k, n) / rho_;
                            }
                        }
                        // Forcing to wave profile instead
                        if (has_outprofile) {
                            // Same steps as in wave generation region
                            amrex::Real new_vof =
                                (1. - Gamma) * target_volfrac(i, j, k) +
                                Gamma * volfrac(i, j, k);
                            new_vof =
                                (new_vof > 1. - vof_tiny)
                                    ? 1.0
                                    : (new_vof < vof_tiny ? 0.0 : new_vof);
                            const amrex::Real dvf = new_vof - volfrac(i, j, k);
                            volfrac(i, j, k) += dvf;
                            const amrex::Real fvel_liq =
                                (target_volfrac(i, j, k) > vof_tiny) ? 1.0
                                                                     : 0.0;
                            amrex::Real rho_ = rho1 * volfrac(i, j, k) +
                                               rho2 * (1.0 - volfrac(i, j, k));
                            for (int n = 0; n < vel.ncomp; ++n) {
                                amrex::Real vel_liq = vel(i, j, k, n);
                                const amrex::Real dvel_liq =
                                    ((1. - Gamma) * target_vel(i, j, k, n) +
                                     Gamma * vel_liq) -
                                    vel_liq;
                                vel_liq += rampf * fvel_liq * dvel_liq;
                                amrex::Real integrated_vel_liq =
                                    volfrac(i, j, k) * vel_liq;
                                integrated_vel_liq +=
                                    rampf * fvel_liq * amrex::max(0.0, dvf) *
                                    (target_vel(i, j, k, n) - vel(i, j, k, n));
                                vel(i, j, k, n) =
                                    (rho1 * integrated_vel_liq +
                                     rho2 * (1. - volfrac(i, j, k)) *
                                         vel(i, j, k, n)) /
                                    rho_;
                            }
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
    vof.fillpatch(time);
    velocity.fillpatch(time);
    density.fillpatch(time);
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
