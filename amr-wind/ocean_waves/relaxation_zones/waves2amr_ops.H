#ifndef W2A_WAVES_OPS_H
#define W2A_WAVES_OPS_H

#include "amr-wind/ocean_waves/relaxation_zones/W2AWaves.H"
#include "amr-wind/ocean_waves/OceanWavesOps.H"
#include "amr-wind/ocean_waves/OceanWaves.H"
#include "amr-wind/ocean_waves/relaxation_zones/relaxation_zones_ops.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/equation_systems/BCOps.H"
#include "AMReX_MultiFabUtil.H"

namespace amr_wind::ocean_waves::ops {

template <>
struct ReadInputsOp<W2AWaves>
{
    void operator()(
        W2AWaves::DataType& data, const ::amr_wind::utils::MultiParser& pp)
    {
// Check for W2A initialization
#ifndef AMR_WIND_USE_W2A
        // Assert Waves2AMR must be used for initial condition file
        amrex::Abort(
            "ocean_waves/W2AWaves: AMR-Wind was not built with Waves2AMR "
            "support; associated wave data cannot be processed for relaxation "
            "zones.");
#endif

        auto& wdata = data.meta();
        // cppcheck-suppress constVariableReference
        auto& info = data.info();
        relaxation_zones::read_inputs(wdata, info, pp);

        pp.get("HOS_modes_filename", wdata.modes_file);
        pp.query("HOS_init_timestep", wdata.ntime);
        if (!pp.contains("HOS_init_timestep")) {
            pp.query("HOS_init_time", wdata.t_winit);
        }

#ifdef AMR_WIND_USE_W2A

        amrex::Vector<amrex::Real> prob_lo_input(AMREX_SPACEDIM);
        amrex::ParmParse pp_geom("geometry");
        pp_geom.getarr("prob_lo", prob_lo_input);

        // Read user inputs to height vector
        int nheights = 0;
        amrex::Real dz0 = 0.;
        pp.get("number_interp_points_in_z", nheights);
        pp.get("interp_spacing_at_surface", dz0);

        // Initialize mode reader
        wdata.rmodes.initialize(wdata.modes_file);

        // Get dt of HOS data
        wdata.dt_modes = wdata.rmodes.get_dtout();

        // Get initial time and timestep synced
        if (wdata.t_winit > 0.0) {
            // If initial time was specified
            // Get time index near requested time
            wdata.ntime = wdata.rmodes.time2step(wdata.t_winit, wdata.ntime);
            // Sync time to time index
            wdata.t_winit = wdata.dt_modes * wdata.ntime;
        } else {
            // If initial timestep is being used
            wdata.t_winit = wdata.dt_modes * wdata.ntime;
        }

        // Initialize variables to store modes
        int vsize = wdata.rmodes.get_vector_size();
        double initval = 0.0;
        wdata.mX.resize(vsize, initval);
        wdata.mY.resize(vsize, initval);
        wdata.mZ.resize(vsize, initval);
        wdata.mFS.resize(vsize, initval);

        // Get dimensions of data
        wdata.n0 = wdata.rmodes.get_first_dimension();
        wdata.n1 = wdata.rmodes.get_second_dimension();
        // Get resolution
        wdata.dx0 = wdata.rmodes.get_xlen() / wdata.n0;
        wdata.dx1 = wdata.rmodes.get_ylen() / wdata.n1;
        // Get depth
        const amrex::Real depth = wdata.rmodes.get_depth();
        // Get dimensional length
        wdata.dimL = wdata.rmodes.get_L();

        // Warning if depth does not correspond to simulation
        if (std::abs(depth - (wdata.zsl - prob_lo_input[2])) > 1e-3 * depth) {
            amrex::Print()
                << "WARNING: Mismatch between water depths from AMR-Wind "
                   "domain and HOS data interpreted by Waves2AMR";
        }

        // Allocate pointers for FFTW
        wdata.eta_mptr = modes_hosgrid::allocate_complex(wdata.n0, wdata.n1);
        wdata.u_mptr = modes_hosgrid::allocate_complex(wdata.n0, wdata.n1);
        wdata.v_mptr = modes_hosgrid::allocate_complex(wdata.n0, wdata.n1);
        wdata.w_mptr = modes_hosgrid::allocate_complex(wdata.n0, wdata.n1);

        // Set up plan for FFTW
        wdata.plan =
            modes_hosgrid::plan_ifftw(wdata.n0, wdata.n1, wdata.eta_mptr);

        // Create height vector for velocity mode conversion before
        // interpolation, with prob_lo as bottom
        int flag = interp_to_mfab::create_height_vector(
            wdata.hvec, nheights, dz0, wdata.zsl, prob_lo_input[2]);
        // Fail if flag indicates it should
        if (flag > 0) {
            amrex::Abort(
                "Waves2AMR ReadInputsOp: create_height_vector error, failure "
                "code " +
                std::to_string(flag));
        }

        // If init_wave_field is activated and initialization will be done, get
        // modes on every processor
        if (wdata.init_wave_field && data.sim().time().time_index() == 0) {
            wdata.rmodes.get_data(
                wdata.ntime, wdata.mX, wdata.mY, wdata.mZ, wdata.mFS);

            // Convert modes to spatial data
            modes_hosgrid::copy_complex(
                wdata.n0, wdata.n1, wdata.mFS, wdata.eta_mptr);
            wdata.sp_eta_vec.resize(wdata.n0 * wdata.n1, 0.0);
            modes_hosgrid::populate_hos_eta(
                wdata.rmodes, wdata.plan, wdata.eta_mptr, wdata.sp_eta_vec);
            // Mesh is not yet created, so get data at every height
            int n_hts = wdata.hvec.size();
            wdata.sp_u_vec.resize(wdata.n0 * wdata.n1 * n_hts);
            wdata.sp_v_vec.resize(wdata.n0 * wdata.n1 * n_hts);
            wdata.sp_w_vec.resize(wdata.n0 * wdata.n1 * n_hts);
            for (int iht = 0; iht < n_hts; ++iht) {
                // Get sample height
                amrex::Real ht = wdata.hvec[iht];
                // Sample velocity
                modes_hosgrid::populate_hos_vel(
                    wdata.rmodes, ht, wdata.mX, wdata.mY, wdata.mZ, wdata.plan,
                    wdata.u_mptr, wdata.v_mptr, wdata.w_mptr, wdata.sp_u_vec,
                    wdata.sp_v_vec, wdata.sp_w_vec, iht * wdata.n0 * wdata.n1);
            }
        }

#endif

        // Declare fields for HOS
        auto& w2a_levelset =
            data.sim().repo().declare_field("w2a_levelset", 1, 3, 1);
        // cppcheck-suppress constVariableReference
        auto& w2a_velocity = data.sim().repo().declare_field(
            "w2a_velocity", AMREX_SPACEDIM, 3, 1);
        // cppcheck-suppress constVariableReference

        // Extrapolation can work well when finer data is available
        w2a_levelset.set_default_fillpatch_bc(data.sim().time());
        w2a_velocity.set_default_fillpatch_bc(data.sim().time());
    }
}; // namespace ops

template <>
struct InitDataOp<W2AWaves>
{
    void
    operator()(W2AWaves::DataType& data, int level, const amrex::Geometry& geom)
    {
        auto& wdata = data.meta();
        bool init_waves = wdata.init_wave_field;

        auto& sim = data.sim();
        // cppcheck-suppress constVariable

        auto& levelset = sim.repo().get_field("levelset");
        auto& velocity = sim.repo().get_field("velocity");

        const auto& problo = geom.ProbLoArray();
        const auto& dx = geom.CellSizeArray();

#ifdef AMR_WIND_USE_W2A

        // Blank initialization if asked for
        if (!init_waves) {
            // Loop to populate field data
            for (amrex::MFIter mfi(levelset(level)); mfi.isValid(); ++mfi) {
                auto phi = levelset(level).array(mfi);
                auto vel = velocity(level).array(mfi);
                const auto& vbx = mfi.validbox();
                const auto& gbx = grow(vbx, 3);
                const amrex::Real zsl = wdata.zsl;
                amrex::ParallelFor(
                    gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                        phi(i, j, k) = zsl - z;
                        vel(i, j, k, 0) = 0.0;
                        vel(i, j, k, 1) = 0.0;
                        vel(i, j, k, 2) = 0.0;
                    });
            }
            return;
        }

        // Set t_last to 0.0 to signify information read in
        wdata.t_last = 0.0;

        // indvec is complete upon initialization (all heights every proc)
        amrex::Vector<int> indvec;
        indvec.resize(wdata.hvec.size());
        for (int n = 0; n < indvec.size(); ++n) {
            indvec[n] = n;
        }
        // Interpolate to MultiFabs (one level at a time)
        interp_to_mfab::interp_eta_to_levelset_multifab(
            wdata.n0, wdata.n1, wdata.dx0, wdata.dx1, wdata.zsl,
            wdata.sp_eta_vec, levelset(level), problo, dx);
        interp_to_mfab::interp_velocity_to_multifab(
            wdata.n0, wdata.n1, wdata.dx0, wdata.dx1, indvec, wdata.hvec,
            wdata.sp_u_vec, wdata.sp_v_vec, wdata.sp_w_vec, velocity(level),
            problo, dx);

        // Copy to ow fields for future interpolation
        auto& ow_levelset = sim.repo().get_field("ow_levelset");
        auto& ow_velocity = sim.repo().get_field("ow_velocity");
        amrex::MultiFab::Copy(ow_levelset(level), levelset(level), 0, 0, 1, 3);
        amrex::MultiFab::Copy(
            ow_velocity(level), velocity(level), 0, 0, AMREX_SPACEDIM, 3);

#else
        amrex::ignore_unused(init_waves, problo, dx);
#endif
    }
}; // namespace ocean_waves

template <>
struct UpdateRelaxZonesOp<W2AWaves>
{
    void operator()(W2AWaves::DataType& data)
    {
        auto& wdata = data.meta();

        auto& sim = data.sim();

        // cppcheck-suppress constVariable
        // Nudge the solution toward where it should be
        amrex::Real time = sim.time().new_time();

        // Update ow fields every time
        auto& m_ow_levelset = sim.repo().get_field("ow_levelset");
        // cppcheck-suppress constVariable
        auto& m_ow_velocity = sim.repo().get_field("ow_velocity");
        // Update HOS fields when necessary
        auto& w2a_levelset = sim.repo().get_field("w2a_levelset");
        auto& w2a_velocity = sim.repo().get_field("w2a_velocity");

        auto nlevels = sim.repo().num_active_levels();
        auto geom = sim.mesh().Geom();

        // Get value for time interpolation
        amrex::Real t_last = wdata.t_last;

#ifdef AMR_WIND_USE_W2A
        // Check if new HOS data needs to be read
        bool read_flag = false;
        // Check if time indicates reading must take place
        int new_ntime =
            wdata.rmodes.time2step(time + wdata.t_winit, wdata.ntime);
        if (new_ntime != wdata.ntime) {
            // New data is needed - reading should happen
            read_flag = true;
            // Time index for reading
            wdata.ntime = new_ntime;
            // Sim time to go with recorded data
            wdata.t = new_ntime * wdata.dt_modes + wdata.t_winit;
        }
        // Check if reading must take place for other reasons
        if (wdata.t_last < -1e-10) {
            // Signifies initialization from scratch without waves or a restart
            read_flag = true;
            // Resizing needs to happen for the first time
            wdata.resize_flag = true;

            if (std::abs(sim.time().current_time()) < 1e-10) {
                // Start from scratch
                // NEED: how do I start from scratch without the original data?
            } else {
                // Restart
                // NEED: how do I restart properly? need old ow fields?
            }

        } else if (std::abs(wdata.t_last) < 1e-10) {
            // Signifies initialization with waves
            read_flag = true;
            // Resizing needs to happen for the first time
            wdata.resize_flag = true;

            // levelset and velocity fields are up-to-date at t=0
            // interpolation is ready to go
        }
        // Record latest time as 'last' for next timestep
        wdata.t_last = time;

        // Resizing (assuming reading is taking place) must happen after regrid
        if (wdata.regrid_occurred) {
            // resize_flag remains true until resizing occurs, but
            // regrid_occurred resets every timestep
            wdata.resize_flag = true;
        }

        // Read HOS data if necessary based on time
        if (read_flag) {

            if (wdata.resize_flag) {
                // Reset flag
                wdata.resize_flag = false;
                // Get heights for this processor
                int flag = interp_to_mfab::get_local_height_indices(
                    wdata.indvec, wdata.hvec, m_ow_velocity.vec_ptrs(), geom);

                // Check lateral bounds (in x)
                const int dir = 0;
                flag += 1 - interp_to_mfab::check_lateral_overlap_lo(
                                wdata.gen_length, dir, m_ow_velocity.vec_ptrs(),
                                geom);
                if (wdata.has_outprofile) {
                    flag += 1 - interp_to_mfab::check_lateral_overlap_hi(
                                    wdata.beach_length, dir,
                                    m_ow_velocity.vec_ptrs(), geom);
                }

                if (flag > 0) {
                    // No overlapping with spatial data or no overlapping with
                    // relaxation zones, interpolation can be skipped
                    wdata.do_interp = false;
                } else {
                    // Do resizing
                    wdata.sp_eta_vec.resize(wdata.n0 * wdata.n1, 0.0);
                    wdata.sp_u_vec.resize(
                        wdata.n0 * wdata.n1 * wdata.indvec.size());
                    wdata.sp_v_vec.resize(
                        wdata.n0 * wdata.n1 * wdata.indvec.size());
                    wdata.sp_w_vec.resize(
                        wdata.n0 * wdata.n1 * wdata.indvec.size());
                }
            }

            // Skip work if no interpolation is needed
            if (!wdata.do_interp) {
                return;
            }

            // Get data from modes
            wdata.rmodes.get_data(
                wdata.ntime, wdata.mX, wdata.mY, wdata.mZ, wdata.mFS);

            // Convert to spatial data in vectors
            modes_hosgrid::copy_complex(
                wdata.n0, wdata.n1, wdata.mFS, wdata.eta_mptr);
            modes_hosgrid::populate_hos_eta(
                wdata.rmodes, wdata.plan, wdata.eta_mptr, wdata.sp_eta_vec);

            for (int iht = 0; iht < wdata.indvec.size(); ++iht) {
                // Get sample height
                amrex::Real ht = wdata.hvec[wdata.indvec[iht]];
                // Sample velocity
                modes_hosgrid::populate_hos_vel(
                    wdata.rmodes, ht, wdata.mX, wdata.mY, wdata.mZ, wdata.plan,
                    wdata.u_mptr, wdata.v_mptr, wdata.w_mptr, wdata.sp_u_vec,
                    wdata.sp_v_vec, wdata.sp_w_vec, iht * wdata.n0 * wdata.n1);
            }

            // Interpolate to fields (vector of MultiFabs)
            interp_to_mfab::interp_eta_to_levelset_field(
                wdata.n0, wdata.n1, wdata.dx0, wdata.dx1, wdata.zsl,
                wdata.sp_eta_vec, w2a_levelset.vec_ptrs(), geom);
            interp_to_mfab::interp_velocity_to_field(
                wdata.n0, wdata.n1, wdata.dx0, wdata.dx1, wdata.indvec,
                wdata.hvec, wdata.sp_u_vec, wdata.sp_v_vec, wdata.sp_w_vec,
                w2a_velocity.vec_ptrs(), geom);

            // Average down to get fine information on coarse grid where
            // possible (may be unnecessary)
            for (int lev = nlevels - 1; lev > 0; --lev) {
                amrex::average_down(
                    w2a_velocity(lev), w2a_velocity(lev - 1), 0, AMREX_SPACEDIM,
                    sim.mesh().refRatio(lev - 1));
                amrex::average_down(
                    w2a_levelset(lev), w2a_levelset(lev - 1), 0, 1,
                    sim.mesh().refRatio(lev - 1));
            }
            // Fill in across internal boundaries
            w2a_velocity.fillpatch(0.0);
            w2a_levelset.fillpatch(0.0);
        }
#endif

        // Temporally interpolate at every timestep to get target solution
        for (int lev = 0; lev < nlevels; ++lev) {
            for (amrex::MFIter mfi(m_ow_levelset(lev)); mfi.isValid(); ++mfi) {
                auto phi = m_ow_levelset(lev).array(mfi);
                auto vel = m_ow_velocity(lev).array(mfi);
                auto W2A_phi = w2a_levelset(lev).array(mfi);
                auto W2A_vel = w2a_velocity(lev).array(mfi);

                const amrex::Real W2A_t = wdata.t;
                const auto& gbx = mfi.growntilebox(3);
                amrex::ParallelFor(
                    gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        // Interpolate temporally every time
                        phi(i, j, k) += (W2A_phi(i, j, k) - phi(i, j, k)) *
                                        (time - t_last) /
                                        (W2A_t - t_last + 1e-16);
                        vel(i, j, k, 0) +=
                            (W2A_vel(i, j, k, 0) - vel(i, j, k, 0)) *
                            (time - t_last) / (W2A_t - t_last + 1e-16);
                        vel(i, j, k, 1) +=
                            (W2A_vel(i, j, k, 1) - vel(i, j, k, 1)) *
                            (time - t_last) / (W2A_t - t_last + 1e-16);
                        vel(i, j, k, 2) +=
                            (W2A_vel(i, j, k, 2) - vel(i, j, k, 2)) *
                            (time - t_last) / (W2A_t - t_last + 1e-16);
                    });
            }
        }
    }
};

} // namespace amr_wind::ocean_waves::ops

#endif /* WAVES2AMR_OPS_H */
