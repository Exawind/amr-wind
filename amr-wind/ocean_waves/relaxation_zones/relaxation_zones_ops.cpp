#include "amr-wind/ocean_waves/relaxation_zones/relaxation_zones_ops.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "amr-wind/core/MultiParser.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/io_utils.H"

#include "amr-wind/fvm/gradient.H"
#include "amr-wind/core/field_ops.H"

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

    pp.query("current", wdata.current);

    wdata.has_ramp = pp.contains("timeramp_period");
    if (wdata.has_ramp) {
        pp.get("timeramp_period", wdata.ramp_period);
    }

    amrex::ParmParse pp_multiphase("MultiPhase");
    pp_multiphase.add("water_level", wdata.zsl);
}

void init_data_structures(RelaxZonesBaseData& /*unused*/) {}

void update_target_vof(CFDSim& sim)
{
    const int nlevels = sim.repo().num_active_levels();
    const auto& ow_levelset = sim.repo().get_field("ow_levelset");
    auto& ow_vof = sim.repo().get_field("ow_vof");
    const auto& geom = sim.mesh().Geom();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = geom[lev].CellSizeArray();
        const auto target_phi = ow_levelset(lev).const_arrays();
        auto target_volfrac = ow_vof(lev).arrays();
        const amrex::Real eps = 2. * std::cbrt(dx[0] * dx[1] * dx[2]);
        amrex::ParallelFor(
            ow_vof(lev), amrex::IntVect(2),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                target_volfrac[nbx](i, j, k) =
                    multiphase::levelset_to_vof(i, j, k, eps, target_phi[nbx]);
            });
    }
    amrex::Gpu::synchronize();
}

void apply_relaxation_zones(CFDSim& sim, const RelaxZonesBaseData& wdata)
{
    const int nlevels = sim.repo().num_active_levels();
    const auto& ow_vof = sim.repo().get_field("ow_vof");
    const auto& ow_vel = sim.repo().get_field("ow_velocity");
    const auto& geom = sim.mesh().Geom();

    const auto& mphase = sim.physics_manager().get<MultiPhase>();
    const amrex::Real rho1 = mphase.rho1();
    const amrex::Real rho2 = mphase.rho2();
    constexpr amrex::Real vof_tiny = 1e-12;

    // Get time
    const auto& time = sim.time().new_time();
    const amrex::Real rampf =
        (wdata.has_ramp) ? utils::ramp(time, wdata.ramp_period) : 1.0;

    auto& vof = sim.repo().get_field("vof");
    auto& velocity = sim.repo().get_field("velocity");
    auto& density = sim.repo().get_field("density");

    amr_wind::IntField* terrain_blank_ptr{nullptr};
    amr_wind::IntField* terrain_drag_ptr{nullptr};
    bool terrain_exists{false};
    // Get fields to prevent forcing in or near underwater terrain
    if (sim.repo().int_field_exists("terrain_blank")) {
        terrain_blank_ptr = &sim.repo().get_int_field("terrain_blank");
        terrain_drag_ptr = &sim.repo().get_int_field("terrain_drag");
        terrain_exists = true;
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();
        const auto& probhi = geom[lev].ProbHiArray();
        auto vel_arrs = velocity(lev).arrays();
        auto rho_arrs = density(lev).arrays();
        auto volfrac_arrs = vof(lev).arrays();
        const auto target_volfrac_arrs = ow_vof(lev).const_arrays();
        const auto target_vel_arrs = ow_vel(lev).const_arrays();

        const auto terrain_blank_flags =
            terrain_exists ? (*terrain_blank_ptr)(lev).const_arrays()
                           : amrex::MultiArray4<int const>();
        const auto terrain_drag_flags =
            terrain_exists ? (*terrain_drag_ptr)(lev).const_arrays()
                           : amrex::MultiArray4<int const>();

        const amrex::Real gen_length = wdata.gen_length;
        const amrex::Real beach_length = wdata.beach_length;
        const amrex::Real beach_length_factor = wdata.beach_length_factor;
        const amrex::Real zsl = wdata.zsl;
        const amrex::Real current = wdata.current;
        const bool has_beach = wdata.has_beach;
        const bool has_outprofile = wdata.has_outprofile;

        amrex::ParallelFor(
            velocity(lev), amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const amrex::Real x = amrex::min(
                    amrex::max(problo[0] + (i + 0.5) * dx[0], problo[0]),
                    probhi[0]);
                const amrex::Real z = amrex::min(
                    amrex::max(problo[2] + (k + 0.5) * dx[2], problo[2]),
                    probhi[2]);

                auto vel = vel_arrs[nbx];
                auto rho = rho_arrs[nbx];
                auto volfrac = volfrac_arrs[nbx];
                const auto target_volfrac = target_volfrac_arrs[nbx];
                const auto target_vel = target_vel_arrs[nbx];

                bool in_or_near_terrain{false};
                if (terrain_exists) {
                    in_or_near_terrain =
                        (terrain_blank_flags[nbx](i, j, k) == 1 ||
                         terrain_drag_flags[nbx](i, j, k) == 1);
                }

                // Generation region
                if (x <= problo[0] + gen_length && !in_or_near_terrain) {
                    const amrex::Real Gamma =
                        utils::gamma_generate(x - problo[0], gen_length);
                    // Get bounded new vof, incorporate with increment
                    amrex::Real new_vof = utils::combine_linear(
                        Gamma, target_volfrac(i, j, k), volfrac(i, j, k));
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
                            utils::combine_linear(
                                Gamma, target_vel(i, j, k, n), vel_liq) -
                            vel_liq;
                        vel_liq += rampf * fvel_liq * dvel_liq;
                        // If liquid was added, that liquid has target_vel
                        amrex::Real integrated_vel_liq =
                            volfrac(i, j, k) * vel_liq;
                        integrated_vel_liq +=
                            rampf * fvel_liq * amrex::max(0.0, dvf) *
                            (target_vel(i, j, k, n) - vel(i, j, k, n));
                        // Update overall velocity using momentum
                        vel(i, j, k, n) =
                            (rho1 * integrated_vel_liq +
                             rho2 * (1. - volfrac(i, j, k)) * vel(i, j, k, n)) /
                            rho_;
                    }
                }
                // Outlet region
                if (x + beach_length >= probhi[0] && !in_or_near_terrain) {
                    const amrex::Real Gamma = utils::gamma_absorb(
                        x - (probhi[0] - beach_length), beach_length,
                        beach_length_factor);
                    // Numerical beach (sponge layer)
                    if (has_beach) {
                        // Get bounded new vof, save increment
                        amrex::Real new_vof = utils::combine_linear(
                            Gamma, utils::free_surface_to_vof(zsl, z, dx[2]),
                            volfrac(i, j, k));
                        new_vof = (new_vof > 1. - vof_tiny)
                                      ? 1.0
                                      : (new_vof < vof_tiny ? 0.0 : new_vof);
                        const amrex::Real dvf = new_vof - volfrac(i, j, k);
                        volfrac(i, j, k) = new_vof;
                        // Conserve momentum when density changes
                        amrex::Real rho_ = rho1 * volfrac(i, j, k) +
                                           rho2 * (1.0 - volfrac(i, j, k));
                        // Target vel is current, 0, 0
                        amrex::Real out_vel{0.0};
                        for (int n = 0; n < vel.ncomp; ++n) {
                            if (n == 0) {
                                out_vel = current;
                            } else {
                                out_vel = 0.0;
                            }
                            const amrex::Real vel_liq = utils::combine_linear(
                                Gamma, out_vel, vel(i, j, k, n));
                            amrex::Real integrated_vel_liq =
                                volfrac(i, j, k) * vel_liq;
                            integrated_vel_liq += amrex::max(0.0, dvf) *
                                                  (out_vel - vel(i, j, k, n));
                            vel(i, j, k, n) = (rho1 * integrated_vel_liq +
                                               rho2 * (1. - volfrac(i, j, k)) *
                                                   vel(i, j, k, n)) /
                                              rho_;
                        }
                    }
                    // Forcing to wave profile instead
                    if (has_outprofile) {
                        // Same steps as in wave generation region
                        amrex::Real new_vof = utils::combine_linear(
                            Gamma, target_volfrac(i, j, k), volfrac(i, j, k));
                        new_vof = (new_vof > 1. - vof_tiny)
                                      ? 1.0
                                      : (new_vof < vof_tiny ? 0.0 : new_vof);
                        const amrex::Real dvf = new_vof - volfrac(i, j, k);
                        volfrac(i, j, k) += dvf;
                        const amrex::Real fvel_liq =
                            (target_volfrac(i, j, k) > vof_tiny) ? 1.0 : 0.0;
                        amrex::Real rho_ = rho1 * volfrac(i, j, k) +
                                           rho2 * (1.0 - volfrac(i, j, k));
                        for (int n = 0; n < vel.ncomp; ++n) {
                            amrex::Real vel_liq = vel(i, j, k, n);
                            const amrex::Real dvel_liq =
                                utils::combine_linear(
                                    Gamma, target_vel(i, j, k, n), vel_liq) -
                                vel_liq;
                            vel_liq += rampf * fvel_liq * dvel_liq;
                            amrex::Real integrated_vel_liq =
                                volfrac(i, j, k) * vel_liq;
                            integrated_vel_liq +=
                                rampf * fvel_liq * amrex::max(0.0, dvf) *
                                (target_vel(i, j, k, n) - vel(i, j, k, n));
                            vel(i, j, k, n) = (rho1 * integrated_vel_liq +
                                               rho2 * (1. - volfrac(i, j, k)) *
                                                   vel(i, j, k, n)) /
                                              rho_;
                        }
                    }
                }

                // Make sure that density is updated before entering the
                // solution
                rho(i, j, k) =
                    rho1 * volfrac(i, j, k) + rho2 * (1. - volfrac(i, j, k));
            });
    }
    amrex::Gpu::synchronize();

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
