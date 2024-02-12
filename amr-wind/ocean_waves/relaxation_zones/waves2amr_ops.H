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

#ifdef AMR_WIND_USE_W2A
namespace {
int evaluate_read_resize(
    int& ntime,
    bool& read_flag,
    bool& resize_flag,
    amrex::Real& wtime,
    amrex::Real& t_last,
    const int new_ntime,
    const amrex::Real wtinit,
    const amrex::Real wdt,
    const amrex::Real newtime)
{
    // Flag to indicate that data must be read twice for interpolation
    int double_data = 0;
    // Check if time indicates reading must take place
    if (new_ntime != ntime) {
        // New data is needed - reading should happen
        read_flag = true;
        // If time index has proceeded more than 1 step
        if (new_ntime > ntime + 1) {
            // Double reading is necessary
            double_data = 1;
        }
        // Time index for reading
        ntime = new_ntime;
        // Sim time to go with recorded data
        wtime = new_ntime * wdt + wtinit;
        // If double reading is deemed necessary, check for convenience
        if (double_data == 1 && std::abs(wtime - newtime) <= 1e-10) {
            // Reading can be done just once, w2a fields replace ow fields
            double_data = 2;
        }
    }
    // Check if reading must take place for other reasons
    if (t_last < -1e-10) {
        // Signifies initialization from scratch without waves or a restart
        read_flag = true;
        // Resizing needs to happen for the first time
        resize_flag = true;

        // Confirm that new time is not coincident with modes time step
        if (std::abs(wtime - newtime) > 1e-10) {
            // Data must be read before and after wtime
            double_data = 1;
        } else {
            // Data need only be read once, w2a fields will replace ow fields
            double_data = 2;
        }

    } else if (std::abs(t_last) < 1e-10) {
        // Signifies initialization with waves
        read_flag = true;
        // Resizing needs to happen for the first time
        resize_flag = true;

        // levelset and velocity fields are up-to-date at t=0
        // interpolation is ready to go
    }
    // Record latest time as 'last' for next timestep
    t_last = newtime;
    // Return flag regarding double reading
    return double_data;
}

void postprocess_velocity_mfab_liquid(
    amrex::MultiFab& vel_mfab,
    amrex::MultiFab& lvs_mfab,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
{
    for (amrex::MFIter mfi(vel_mfab); mfi.isValid(); ++mfi) {
        auto vel = vel_mfab.array(mfi);
        auto phi = lvs_mfab.const_array(mfi);
        const auto& vbx = mfi.validbox();
        const auto& gbx = grow(vbx, 3);
        amrex::ParallelFor(
            gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Set velocity to zero if no liquid present
                if (phi(i, j, k) + 0.5 * dx[2] < 0.0) {
                    vel(i, j, k, 0) = 0.0;
                    vel(i, j, k, 1) = 0.0;
                    vel(i, j, k, 2) = 0.0;
                }
            });
    }
}

void postprocess_velocity_field_liquid(
    amr_wind::Field& vel_field,
    amr_wind::Field& lvs_field,
    amrex::Vector<amrex::Geometry> geom_all)
{
    int nlevels = vel_field.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx_lev = geom_all[lev].CellSizeArray();
        postprocess_velocity_mfab_liquid(
            vel_field(lev), lvs_field(lev), dx_lev);
    }
}

int update_offset_timestep(const int ntime, const int n0)
{
    // Offending timestep (goes too far): ntime + offset
    // Farthest back timesteps are permitted to go: n0
    // Subtract offset by offending timestep, add back lower limit
    // new offset = offset - (ntime + offset) + n0
    return (-ntime + n0);
}

void populate_fields_all_levels(
    amr_wind::ocean_waves::W2AWaves::MetaType& wdata,
    amrex::Vector<amrex::Geometry> geom_all,
    amr_wind::Field& lvs_field,
    amr_wind::Field& vel_field,
    int ntime_off = 0)
{

    // Get data from modes
    bool no_EOF = wdata.rmodes.get_data(
        wdata.ntime + wdata.n_offset + ntime_off, wdata.mX, wdata.mY, wdata.mZ,
        wdata.mFS);
    // Navigate when end of file is reached
    if (!no_EOF) {
        // End of file detected, reset reading
        wdata.n_offset = update_offset_timestep(wdata.ntime, wdata.n_winit);
        // Print warning to screen
        amrex::Print() << "WARNING (waves2amr_ops): end of mode data file "
                          "detected, resetting to beginning of mode data.\n";
        // Read data again, now from a valid timestep
        no_EOF = wdata.rmodes.get_data(
            wdata.ntime + wdata.n_offset + ntime_off, wdata.mX, wdata.mY,
            wdata.mZ, wdata.mFS);
        // If no valid data is detected at this point, abort
        if (!no_EOF) {
            amrex::Abort(
                "waves2amr_ops: end of mode data file detected after "
                "resetting to beginning; please evaluate HOS_init_time "
                "or HOS_init_timestep and check the length of the mode "
                "file.");
        }
    }

    // Convert to spatial data in vectors
    modes_hosgrid::copy_complex(wdata.n0, wdata.n1, wdata.mFS, wdata.eta_mptr);
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
        wdata.n0, wdata.n1, wdata.dx0, wdata.dx1, wdata.zsl, wdata.sp_eta_vec,
        lvs_field.vec_ptrs(), geom_all);
    interp_to_mfab::interp_velocity_to_field(
        wdata.n0, wdata.n1, wdata.dx0, wdata.dx1, wdata.indvec, wdata.hvec,
        wdata.sp_u_vec, wdata.sp_v_vec, wdata.sp_w_vec, vel_field.vec_ptrs(),
        geom_all);

    // Zero velocity in pure gas cells
    postprocess_velocity_field_liquid(vel_field, lvs_field, geom_all);
}

} // namespace
#endif

namespace amr_wind::ocean_waves::ops {

template <>
struct ReadInputsOp<W2AWaves>
{
    void operator()(
        // cppcheck-suppress constParameterReference
        W2AWaves::DataType& data,
        const ::amr_wind::utils::MultiParser& pp)
    {
// Check for W2A initialization
#ifndef AMR_WIND_USE_W2A
        // Assert Waves2AMR must be used for initial condition file
        amrex::Abort(
            "ocean_waves/W2AWaves: AMR-Wind was not built with Waves2AMR "
            "support; associated wave data cannot be processed for relaxation "
            "zones.");

        amrex::ignore_unused(data, pp);
#else
        auto& wdata = data.meta();
        auto& info = data.info();
        relaxation_zones::read_inputs(wdata, info, pp);

        pp.get("HOS_modes_filename", wdata.modes_file);
        pp.query("HOS_init_timestep", wdata.ntime);
        if (!pp.contains("HOS_init_timestep")) {
            pp.query("HOS_init_time", wdata.t_winit);
        }

        amrex::Vector<amrex::Real> prob_lo_input(AMREX_SPACEDIM);
        amrex::ParmParse pp_geom("geometry");
        pp_geom.getarr("prob_lo", prob_lo_input);

        // Read user inputs to height vector
        int nheights = 0;
        int nh_above = 1;
        amrex::Real dz0 = 0.;
        pp.get("number_interp_points_in_z", nheights);
        pp.get("interp_spacing_at_surface", dz0);
        pp.query("number_interp_above_surface", nh_above);

        // Initialize mode reader
        bool file_exists = wdata.rmodes.initialize(wdata.modes_file);

        // Abort if file does not exist
        if (!file_exists) {
            amrex::Abort(
                "Waves2AMR ReadInputsOp: modes file requested does not exist");
        }

        // Get dt of HOS data
        wdata.dt_modes = wdata.rmodes.get_dtout();

        // Get initial time and timestep synced
        if (wdata.t_winit > 0.0) {
            // If initial time was specified
            // Get time index near requested time
            wdata.ntime = wdata.rmodes.time2step(wdata.t_winit, wdata.ntime);
            // Sync time to time index
            wdata.t_winit = wdata.dt_modes * wdata.ntime;
            // Save first timestep
            wdata.n_winit = wdata.ntime;
        } else {
            // If initial timestep is being used
            wdata.t_winit = wdata.dt_modes * wdata.ntime;
            // Save first timestep
            wdata.n_winit = wdata.ntime;
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
        // Get nominal last timestep of data
        wdata.n_wstop =
            (int)((wdata.rmodes.get_Tstop() + 1e-8) / wdata.dt_modes);

        // Check if stop time is exceeded, introduce offset to ntime
        if (wdata.ntime + wdata.n_offset > wdata.n_wstop) {
            // If exceeding stop step, calculate new offset
            wdata.n_offset = update_offset_timestep(wdata.ntime, wdata.n_winit);
            // Print warning to screen
            amrex::Print()
                << "WARNING (waves2amr_ops): available mode data exceeded, "
                   "resetting to beginning of mode data.\n";
        }

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
            wdata.hvec, nheights, dz0, wdata.zsl, prob_lo_input[2], nh_above);
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
            bool no_EOF = wdata.rmodes.get_data(
                wdata.ntime + wdata.n_offset, wdata.mX, wdata.mY, wdata.mZ,
                wdata.mFS);
            if (!no_EOF) {
                // End of file detected, reset reading
                wdata.n_offset =
                    update_offset_timestep(wdata.ntime, wdata.n_winit);
                // Print warning to screen
                amrex::Print()
                    << "WARNING (waves2amr_ops): end of mode data file "
                       "detected, resetting to beginning of mode data.\n";
                // Read data again, now from a valid timestep
                no_EOF = wdata.rmodes.get_data(
                    wdata.ntime + wdata.n_offset, wdata.mX, wdata.mY, wdata.mZ,
                    wdata.mFS);
                // If no valid data is detected at this point, abort
                if (!no_EOF) {
                    amrex::Abort(
                        "waves2amr_ops: end of mode data file detected after "
                        "resetting to beginning; please evaluate HOS_init_time "
                        "or HOS_init_timestep and check the length of the mode "
                        "file.");
                }
            }

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

        // Declare fields for HOS
        auto& w2a_levelset =
            data.sim().repo().declare_field("w2a_levelset", 1, 3, 1);
        auto& w2a_velocity = data.sim().repo().declare_field(
            "w2a_velocity", AMREX_SPACEDIM, 3, 1);

        // Extrapolation can work well when finer data is available
        w2a_levelset.set_default_fillpatch_bc(data.sim().time());
        w2a_velocity.set_default_fillpatch_bc(data.sim().time());
#endif
    }
}; // namespace ops

template <>
struct InitDataOp<W2AWaves>
{
    void
    // cppcheck-suppress constParameterReference
    operator()(W2AWaves::DataType& data, int level, const amrex::Geometry& geom)
    {

#ifdef AMR_WIND_USE_W2A
        auto& wdata = data.meta();
        bool init_waves = wdata.init_wave_field;

        auto& sim = data.sim();

        auto& levelset = sim.repo().get_field("levelset");
        auto& velocity = sim.repo().get_field("velocity");

        const auto& problo = geom.ProbLoArray();
        const auto& dx = geom.CellSizeArray();

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
        // Zero velocity in pure liquid cells
        postprocess_velocity_mfab_liquid(velocity(level), levelset(level), dx);

        // Copy to ow fields for future interpolation
        auto& ow_levelset = sim.repo().get_field("ow_levelset");
        auto& ow_velocity = sim.repo().get_field("ow_velocity");
        amrex::MultiFab::Copy(ow_levelset(level), levelset(level), 0, 0, 1, 3);
        amrex::MultiFab::Copy(
            ow_velocity(level), velocity(level), 0, 0, AMREX_SPACEDIM, 3);
#else
        amrex::ignore_unused(data, level, geom);
#endif
    }
}; // namespace ocean_waves

template <>
struct UpdateRelaxZonesOp<W2AWaves>
{
    // cppcheck-suppress constParameterReference
    void operator()(W2AWaves::DataType& data)
    {

#ifdef AMR_WIND_USE_W2A
        auto& wdata = data.meta();
        auto& sim = data.sim();

        // Nudge the solution toward where it should be
        const amrex::Real newtime = sim.time().new_time();

        // Update ow fields every time
        auto& m_ow_levelset = sim.repo().get_field("ow_levelset");
        auto& m_ow_velocity = sim.repo().get_field("ow_velocity");
        // Update HOS fields when necessary
        auto& w2a_levelset = sim.repo().get_field("w2a_levelset");
        auto& w2a_velocity = sim.repo().get_field("w2a_velocity");

        auto nlevels = sim.repo().num_active_levels();
        auto geom = sim.mesh().Geom();

        // Get value for time interpolation
        amrex::Real t_last = wdata.t_last;

        // Check if new HOS data needs to be read
        bool read_flag = false;
        // Check if time indicates reading must take place
        int new_ntime =
            wdata.rmodes.time2step(newtime + wdata.t_winit, wdata.ntime);
        int double_data = evaluate_read_resize(
            wdata.ntime, read_flag, wdata.resize_flag, wdata.t, wdata.t_last,
            new_ntime, wdata.t_winit, wdata.dt_modes, newtime);
        // Check if stop time is exceeded, introduce offset to ntime
        if (read_flag) {
            // Need to only check when reading is occurring
            if (wdata.ntime + wdata.n_offset > wdata.n_wstop) {
                // If exceeding stop step, calculate new offset
                wdata.n_offset =
                    update_offset_timestep(wdata.ntime, wdata.n_winit);
                // Print warning to screen
                amrex::Print()
                    << "WARNING (waves2amr_ops): available mode data exceeded, "
                       "resetting to beginning of mode data.\n";
            }
        }
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
                // Flags for indicating overlap, assume none at first
                bool flag_z = false;
                bool flag_xlo = false;
                bool flag_xhi = false;
                // Get heights for this processor, check overlap in z
                flag_z =
                    (interp_to_mfab::get_local_height_indices(
                         wdata.indvec, wdata.hvec, m_ow_velocity.vec_ptrs(),
                         geom) == 1);
                // No overlap from heights definitely means no interp

                // Check lateral bounds (in x)
                const int dir = 0;
                flag_xlo =
                    (interp_to_mfab::check_lateral_overlap_lo(
                         wdata.gen_length, dir, m_ow_velocity.vec_ptrs(),
                         geom) == 1);
                // No overlap with gen region means no interp, unless ...
                if (wdata.has_outprofile) {
                    // ... if overlap exists here, needing interp
                    flag_xhi =
                        (interp_to_mfab::check_lateral_overlap_hi(
                             wdata.beach_length, dir, m_ow_velocity.vec_ptrs(),
                             geom) == 1);
                }

                if (flag_z && (flag_xlo || flag_xhi)) {
                    // Interpolation is needed
                    wdata.do_interp = true;
                    // Do resizing
                    wdata.sp_eta_vec.resize(wdata.n0 * wdata.n1, 0.0);
                    wdata.sp_u_vec.resize(
                        wdata.n0 * wdata.n1 * wdata.indvec.size());
                    wdata.sp_v_vec.resize(
                        wdata.n0 * wdata.n1 * wdata.indvec.size());
                    wdata.sp_w_vec.resize(
                        wdata.n0 * wdata.n1 * wdata.indvec.size());
                    // Sizes will remain constant and need for interpolation
                    // will remain until a regrid occurs
                } else {
                    // No overlapping with spatial data or no overlapping with
                    // relaxation zones, interpolation can be skipped
                    wdata.do_interp = false;
                }
            }

            // Skip work if no interpolation is needed
            if (!wdata.do_interp) {
                return;
            }

            // If double read is required, then copy older wave data to ow_
            // fields and modify interpolation parameters to get things right
            if (double_data == 1) {
                populate_fields_all_levels(
                    wdata, geom, m_ow_levelset, m_ow_velocity, -1);

                // Average down to get fine information on coarse grid where
                // possible (may be unnecessary)
                for (int lev = nlevels - 1; lev > 0; --lev) {
                    amrex::average_down(
                        m_ow_velocity(lev), m_ow_velocity(lev - 1), 0,
                        AMREX_SPACEDIM, sim.mesh().refRatio(lev - 1));
                    amrex::average_down(
                        m_ow_levelset(lev), m_ow_levelset(lev - 1), 0, 1,
                        sim.mesh().refRatio(lev - 1));
                }
                // Fill in across internal boundaries (only goes with avg down)
                m_ow_velocity.fillpatch(0.0);
                m_ow_levelset.fillpatch(0.0);

                // Prior t_last (corresponding to ow fields)
                t_last = (wdata.ntime - 1) * wdata.dt_modes;
            } else if (double_data == 2) {
                // Restarting simulation or taking a big step, new time at ntime
                // Initialize ow fields to 0 for time interp, will be replaced
                m_ow_levelset.setVal(0.0);
                m_ow_velocity.setVal(0.0);
                // No modification needed for t_last, leads to interp factor = 1
            }

            // After possible prior read, now read data for this ntime
            populate_fields_all_levels(wdata, geom, w2a_levelset, w2a_velocity);

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
            // Fill in across internal boundaries (only goes with avg down)
            w2a_velocity.fillpatch(0.0);
            w2a_levelset.fillpatch(0.0);
        }

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
                                        (newtime - t_last) /
                                        (W2A_t - t_last + 1e-16);
                        vel(i, j, k, 0) +=
                            (W2A_vel(i, j, k, 0) - vel(i, j, k, 0)) *
                            (newtime - t_last) / (W2A_t - t_last + 1e-16);
                        vel(i, j, k, 1) +=
                            (W2A_vel(i, j, k, 1) - vel(i, j, k, 1)) *
                            (newtime - t_last) / (W2A_t - t_last + 1e-16);
                        vel(i, j, k, 2) +=
                            (W2A_vel(i, j, k, 2) - vel(i, j, k, 2)) *
                            (newtime - t_last) / (W2A_t - t_last + 1e-16);
                    });
            }
        }
#else
        amrex::ignore_unused(data);
#endif
    }
};

} // namespace amr_wind::ocean_waves::ops

#endif /* WAVES2AMR_OPS_H */
