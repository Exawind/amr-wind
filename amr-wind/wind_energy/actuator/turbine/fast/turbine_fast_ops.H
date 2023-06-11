#ifndef TURBINE_FAST_OPS_H
#define TURBINE_FAST_OPS_H

#include "amr-wind/wind_energy/actuator/turbine/fast/TurbineFast.H"
#include "amr-wind/wind_energy/actuator/turbine/turbine_utils.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/wind_energy/actuator/actuator_utils.H"
#include "amr-wind/wind_energy/actuator/FLLCOp.H"

namespace amr_wind {
namespace actuator {
namespace ops {

template <typename SrcTrait>
struct ReadInputsOp<TurbineFast, SrcTrait>
{
    void operator()(TurbineFast::DataType& data, const utils::ActParser& pp)
    {
        // Data common to any turbine actuator simulation
        utils::read_inputs(data.meta(), data.info(), pp);

        auto& tdata = data.meta();

        // Get density value for normalization
        pp.query("density", tdata.density);

        // Initialize OpenFAST specific data
        const auto& tinfo = data.info();
        auto& tf = data.meta().fast_data;
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
            tf.base_pos[i] = static_cast<float>(tinfo.base_pos[i]);

        tf.tlabel = tinfo.label;
        tf.tid_global = tinfo.id;
        tf.num_blades = tdata.num_blades;
        tf.num_pts_blade = tdata.num_pts_blade;
        tf.num_pts_tower = tdata.num_pts_tower;
        tf.dt_cfd = data.sim().time().deltaT();

        pp.get("openfast_start_time", tf.start_time);
        pp.get("openfast_stop_time", tf.stop_time);

        std::string sim_mode = (tf.start_time > 0.0) ? "restart" : "init";
        pp.query("openfast_sim_mode", sim_mode);

        if (sim_mode == "init") {
            tf.sim_mode = ::exw_fast::SimMode::init;
            amrex::Print() << "Initializing turbine:" << tf.tlabel << std::endl;
        } else if (sim_mode == "replay") {
            tf.sim_mode = ::exw_fast::SimMode::replay;
            amrex::Print() << "Replaying turbine:" << tf.tlabel << std::endl;
        } else if (sim_mode == "restart") {
            tf.sim_mode = ::exw_fast::SimMode::restart;
            amrex::Print() << "Restarting turbine:" << tf.tlabel << std::endl;
        } else {
            amrex::Abort(
                "Actuator: Invalid OpenFAST simulation mode: " + sim_mode);
        }

        // If we are using OpenFAST restart file, require that the user provide
        // the path to the checkpoint file.
        if (tf.sim_mode == ::exw_fast::SimMode::restart) {
            pp.get("openfast_restart_file", tf.checkpoint_file);
        } else {
            pp.get("openfast_input_file", tf.input_file);
        }

        const auto& time = data.sim().time();
        tf.chkpt_interval = time.chkpt_interval();

        perform_checks(data);
    }

    void perform_checks(typename TurbineFast::DataType& data)
    {
        const auto& time = data.sim().time();
        // Ensure that we are using fixed timestepping scheme
        AMREX_ALWAYS_ASSERT(!time.adaptive_timestep());
    }
};

template <>
inline void
determine_influenced_procs<TurbineFast>(typename TurbineFast::DataType& data)
{
    auto& info = data.info();
    info.procs =
        utils::determine_influenced_procs(data.sim().mesh(), info.bound_box);

    AMREX_ALWAYS_ASSERT(info.root_proc > -1);
    // During regrid, the influenced processes might have changed and might
    // no longer include the root proc. We insert it back to ensure that it
    // is always present on the list.
    info.procs.insert(info.root_proc);

    const int iproc = amrex::ParallelDescriptor::MyProc();
    auto in_proc = info.procs.find(iproc);
    info.actuator_in_proc = (in_proc != info.procs.end());
    info.sample_vel_in_proc = info.is_root_proc;
}

template <>
inline void determine_root_proc<TurbineFast>(
    typename TurbineFast::DataType& data, amrex::Vector<int>& act_proc_count)
{
    namespace utils = ::amr_wind::actuator::utils;
    auto& info = data.info();
    info.procs =
        utils::determine_influenced_procs(data.sim().mesh(), info.bound_box);

    utils::determine_root_proc(info, act_proc_count);

    // TODO: This function is doing a lot more than advertised by the name.
    // Should figure out a better way to perform the extra work.

    // For OpenFAST we only need velocities sampled in root process
    info.sample_vel_in_proc = info.is_root_proc;

    // Initialize the OpenFAST object and register this turbine in the root
    // process
    if (info.is_root_proc) {
        auto& tdata = data.meta();
        auto& ext_mgr = data.sim().ext_solver_manager();
        ext_mgr.create("OpenFAST", data.sim());
        tdata.fast = &(ext_mgr.get<::exw_fast::FastIface>());
        tdata.fast->register_turbine(tdata.fast_data);
    }
}

template <typename SrcTrait>
struct InitDataOp<TurbineFast, SrcTrait>
{
    void operator()(TurbineFast::DataType& data)
    {
        BL_PROFILE("amr-wind::InitDataOp<TurbineFast>");

        // Ensure that FAST simulation time is set properly before doing any
        // initialization tasks. We perform this check here to account for
        // restart which is only known after reading the checkpoint file.
        check_fast_sim_time(data);

        const auto& info = data.info();
        auto& tdata = data.meta();

        // Initialize our communicator for broadcasting data
        amrex::ParallelDescriptor::Comm_dup(
            amrex::ParallelDescriptor::Communicator(), tdata.tcomm);

        int sz_info[4]{0, 0, 0, 0};
        if (info.is_root_proc) {
            tdata.fast->init_turbine(tdata.fast_data.tid_local);

            const auto& tf = tdata.fast_data;
            sz_info[0] = tf.num_blades;
            sz_info[1] = tf.to_cfd.fx_Len;
            sz_info[2] = tf.from_cfd.u_Len;
            sz_info[3] = tf.num_pts_tower;
        }

        // Broadcast data to everyone
        amrex::ParallelDescriptor::Bcast(
            sz_info, 3, info.root_proc, tdata.tcomm);

        {
            tdata.num_blades = sz_info[0];
            // back calculate what the value per blade for number of points in
            // the openfast data structure
            tdata.num_vel_pts_blade = sz_info[2] / tdata.num_blades - 1;
            data.grid().resize(sz_info[1], sz_info[2]);
            tdata.chord.resize(sz_info[1]);
            tdata.num_pts_tower = sz_info[3];
        }

        tdata.vel_rel.assign(sz_info[1], vs::Vector::zero());

        if (info.is_root_proc) {
            // copy chord data
            int npts = sz_info[1];
            const auto& fchord = tdata.fast_data.to_cfd.forceNodesChord;
            for (int i = 0; i < npts; ++i) {
                tdata.chord[i] = static_cast<amrex::Real>(fchord[i]);
            }
        }

        amrex::ParallelDescriptor::Bcast(
            tdata.chord.data(), tdata.chord.size(), info.root_proc,
            tdata.tcomm);

        make_component_views(data);
        init_epsilon(data);
    }

    void make_component_views(typename TurbineFast::DataType& data)
    {
        auto& grid = data.grid();
        auto& tdata = data.meta();
        const int num_blades = tdata.num_blades;
        const int num_pts_blade = tdata.num_pts_blade;
        const int num_vel_pts_blade = tdata.num_vel_pts_blade;

        for (int ib = 0; ib < num_blades; ++ib) {
            ComponentView cv;

            const auto start = ib * num_pts_blade + 1;
            const auto start_vel = ib * num_vel_pts_blade;
            // clang-format off
            cv.pos = ::amr_wind::utils::slice(
                grid.pos, start, num_pts_blade);
            cv.force = ::amr_wind::utils::slice(
                grid.force, start, num_pts_blade);
            cv.epsilon = ::amr_wind::utils::slice(
                grid.epsilon, start, num_pts_blade);
            cv.orientation = ::amr_wind::utils::slice(
                grid.orientation, start, num_pts_blade);
            cv.chord = ::amr_wind::utils::slice(
                tdata.chord, start, num_pts_blade);
            cv.vel_rel = ::amr_wind::utils::slice(
                tdata.vel_rel, start, num_pts_blade);
            cv.vel= ::amr_wind::utils::slice(
                grid.vel, start_vel, num_vel_pts_blade);
            cv.vel_pos= ::amr_wind::utils::slice(
                grid.vel_pos, start_vel, num_vel_pts_blade);
            // clang-format on

            tdata.blades.emplace_back(std::move(cv));
        }
        if (tdata.num_pts_tower > 0) {
            const int num_pts_tower = tdata.num_pts_tower;
            const int ntwr_start = num_blades * num_pts_blade + 1;
            auto& cv = tdata.tower;

            cv.pos =
                ::amr_wind::utils::slice(grid.pos, ntwr_start, num_pts_tower);
            cv.force =
                ::amr_wind::utils::slice(grid.force, ntwr_start, num_pts_tower);
            cv.epsilon = ::amr_wind::utils::slice(
                grid.epsilon, ntwr_start, num_pts_tower);
            cv.orientation = ::amr_wind::utils::slice(
                grid.orientation, ntwr_start, num_pts_tower);
            cv.chord = ::amr_wind::utils::slice(
                tdata.chord, ntwr_start, num_pts_tower);
        }
        {
            auto& cv = tdata.hub;
            cv.pos = ::amr_wind::utils::slice(grid.pos, 0, 1);
            cv.force = ::amr_wind::utils::slice(grid.force, 0, 1);
            cv.epsilon = ::amr_wind::utils::slice(grid.epsilon, 0, 1);
            cv.orientation = ::amr_wind::utils::slice(grid.orientation, 0, 1);
            cv.chord = ::amr_wind::utils::slice(tdata.chord, 0, 1);
        }
    }

    void init_epsilon(typename TurbineFast::DataType& data)
    {
        auto& tdata = data.meta();

        // Swap order of epsilon based on FAST turbine orientation
        swap_epsilon(tdata.eps_inp);
        swap_epsilon(tdata.eps_min);
        swap_epsilon(tdata.eps_chord);
        swap_epsilon(tdata.eps_tower);

        {
            const auto& cd = tdata.nacelle_cd;
            const auto& area = tdata.nacelle_area;
            const auto eps =
                std::sqrt(2.0 / ::amr_wind::utils::pi() * cd * area);

            auto& nac_eps = data.grid().epsilon[0];
            nac_eps.x() = amrex::max(eps, tdata.eps_min.x());
            nac_eps.y() = amrex::max(eps, tdata.eps_min.y());
            nac_eps.z() = amrex::max(eps, tdata.eps_min.z());
        }

        for (int ib = 0; ib < tdata.num_blades; ++ib) {
            auto& cv = tdata.blades[ib];

            for (int i = 0; i < tdata.num_pts_blade; ++i) {
                const auto eps_crd = tdata.eps_chord * cv.chord[i];

                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    cv.epsilon[i][n] = amrex::max(
                        tdata.eps_min[n], tdata.eps_inp[n], eps_crd[n]);
                }
            }
        }
        {
            auto& cv = tdata.tower;
            for (int i = 0; i < tdata.num_pts_tower; ++i) {
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    cv.epsilon[i][n] = amrex::max(
                        tdata.eps_min[n], tdata.eps_inp[n], tdata.eps_tower[n]);
                }
            }
        }
    }

    inline void swap_epsilon(vs::Vector& eps)
    {
        const auto x = eps.x();
        const auto y = eps.y();
        eps.x() = y;
        eps.y() = x;
    }

    void check_fast_sim_time(typename TurbineFast::DataType& data)
    {
        const auto& time = data.sim().time();

        // Set OpenFAST end time to be at least as long as the CFD time. User
        // can choose a longer duration in input file.
        const amrex::Real stop1 = time.stop_time() > 0.0
                                      ? time.stop_time()
                                      : std::numeric_limits<amrex::Real>::max();
        const amrex::Real stop2 = time.stop_time_index() > -1
                                      ? time.stop_time_index() * time.deltaT()
                                      : std::numeric_limits<amrex::Real>::max();
        const amrex::Real cfd_stop = amrex::min(stop1, stop2);
        const amrex::Real cfd_start = time.current_time();
        const amrex::Real cfd_sim = cfd_stop - cfd_start - 1.0e-6;

        // Ensure that the user specified stop_time is not shorter than CFD sim
        const auto& tf = data.meta().fast_data;
        const amrex::Real fast_sim = tf.stop_time - tf.start_time;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            fast_sim > cfd_sim,
            "OpenFAST simulation time is shorter than AMR-Wind duration");
    }
};

template <typename SrcTrait>
struct UpdatePosOp<TurbineFast, SrcTrait>
{
    void operator()(typename TurbineFast::DataType& data)
    {
        // Return early if this is not the root process for this turbine
        //
        // This is handled in Actuator class, but we add a check here just as a
        // safeguard
        if (!data.info().is_root_proc) return;
        BL_PROFILE("amr-wind::actuator::UpdatePosOp<TurbineFast>");

        const auto& tdata = data.meta();
        const auto& bp = data.info().base_pos;
        const auto& pxvel = tdata.fast_data.to_cfd.pxVel;
        const auto& pyvel = tdata.fast_data.to_cfd.pyVel;
        const auto& pzvel = tdata.fast_data.to_cfd.pzVel;
        auto& vel_pos = data.grid().vel_pos;
        const int nvpts = vel_pos.size();
        for (int i = 0; i < nvpts; ++i) {
            vel_pos[i].x() = static_cast<amrex::Real>(pxvel[i]) + bp.x();
            vel_pos[i].y() = static_cast<amrex::Real>(pyvel[i]) + bp.y();
            vel_pos[i].z() = static_cast<amrex::Real>(pzvel[i]) + bp.z();
        }
    }
};

template <typename SrcTrait>
struct UpdateVelOp<TurbineFast, SrcTrait>
{
    void operator()(typename TurbineFast::DataType& data)
    {
        // Return early if this is not the root process for this turbine
        //
        // This is handled in Actuator class, but we add a check here just as a
        // safeguard
        if (!data.info().is_root_proc) return;
        BL_PROFILE("amr-wind::actuator::UpdateVelOp<TurbineFast>");
        auto& tdata = data.meta();

        auto& from_cfd = tdata.fast_data.from_cfd;
        auto& uvel = from_cfd.u;
        auto& vvel = from_cfd.v;
        auto& wvel = from_cfd.w;
        const auto& vel = data.grid().vel;
        const int nvpts = tdata.vel_rel.size();

        if (!tdata.fllc.empty()) {
            // Compute the relative velocity needed for the FLLC
            for (int i = 0; i < nvpts; ++i) {
                tdata.vel_rel[i][0] =
                    uvel[i] - static_cast<amrex::Real>(
                                  tdata.fast_data.to_cfd.xdotForce[i]);
                tdata.vel_rel[i][1] =
                    vvel[i] - static_cast<amrex::Real>(
                                  tdata.fast_data.to_cfd.ydotForce[i]);
                tdata.vel_rel[i][2] =
                    wvel[i] - static_cast<amrex::Real>(
                                  tdata.fast_data.to_cfd.zdotForce[i]);
            }
            // Loop through each blade and apply the FLLC
            for (int i = 0; i < tdata.num_blades; ++i) {
                FLLCOp()(tdata.blades[i], tdata.fllc[i]);
            }
        }

        for (int i = 0; i < from_cfd.u_Len; ++i) {
            uvel[i] = static_cast<float>(vel[i].x());
            vvel[i] = static_cast<float>(vel[i].y());
            wvel[i] = static_cast<float>(vel[i].z());
        }
    }
};

template <typename SrcTrait>
struct ComputeForceOp<TurbineFast, SrcTrait>
{
    void operator()(typename TurbineFast::DataType& data)
    {
        BL_PROFILE("amr-wind::actuator::ComputeForceOp<TurbineFast>");
        // Advance OpenFAST by specified number of sub-steps
        fast_step(data);
        // Broadcast data to all the processes that contain patches influenced
        // by this turbine
        scatter_data(data);

        const auto& time = data.sim().time();

        auto& tdata = data.meta();
        if (!tdata.fllc.empty()) {
            for (int i = 0; i < tdata.num_blades; ++i) {
                if (!(tdata.fllc[i].initialized) &&
                    (time.current_time() > tdata.fllc[i].fllc_start_time)) {
                    FLLCInit(
                        tdata.fllc[i], tdata.blades[i], tdata.eps_chord[0]);
                }
            }
        }
    }

    void fast_step(typename TurbineFast::DataType& data)
    {
        if (!data.info().is_root_proc) return;

        auto& meta = data.meta();
        auto& tf = data.meta().fast_data;
        if (tf.is_solution0) {
            meta.fast->init_solution(tf.tid_local);
        } else {
            meta.fast->advance_turbine(tf.tid_local);
        }

        // Populate nacelle force into the OpenFAST data structure so that it
        // gets broadcasted to all influenced processes in subsequent scattering
        // of data.
        compute_nacelle_force(data);
    }

    void compute_nacelle_force(typename TurbineFast::DataType& data)
    {
        if (!data.info().is_root_proc) return;

        const auto& cd = data.meta().nacelle_cd;
        const auto& area = data.meta().nacelle_area;
        const auto& cd_area = cd * area;
        const auto& fcfd = data.meta().fast_data.from_cfd;
        const auto& tcfd = data.meta().fast_data.to_cfd;
        const auto& rho = data.meta().density;

        const auto& eps = data.grid().epsilon[0].x();
        vs::Vector vel{fcfd.u[0], fcfd.v[0], fcfd.w[0]};
        amrex::Real correction = 0.0;
        if (eps > 0.0) {
            amrex::Real fac =
                1.0 -
                (cd_area) / (2.0 * ::amr_wind::utils::two_pi() * eps * eps);
            correction = 1.0 / fac;
        }
        amrex::Real coeff =
            0.5 * rho * cd_area * vs::mag(vel) * correction * correction;

        tcfd.fx[0] = static_cast<float>(coeff * fcfd.u[0]);
        tcfd.fy[0] = static_cast<float>(coeff * fcfd.v[0]);
        tcfd.fz[0] = static_cast<float>(coeff * fcfd.w[0]);
    }

    void scatter_data(typename TurbineFast::DataType& data)
    {
        if (!data.info().actuator_in_proc) return;

        // Create an MPI transfer buffer that packs all data in one contiguous
        // array. 3 floats for the position vector, 3 floats for the force
        // vector, and 9 floats for the orientation matrix = 15 floats per
        // actuator node.
        const int dsize = data.grid().pos.size() * 15;
        amrex::Vector<float> buf(dsize);

        // Copy data into MPI send/recv buffer from the OpenFAST data structure.
        // Note, other procs do not have a valid data in those pointers.
        if (data.info().is_root_proc) {
            BL_PROFILE(
                "amr-wind::actuator::ComputeForceOp<TurbineFast>::scatter1");
            const auto& tocfd = data.meta().fast_data.to_cfd;
            auto it = buf.begin();
            std::copy(tocfd.fx, tocfd.fx + tocfd.fx_Len, it);
            std::advance(it, tocfd.fx_Len);
            std::copy(tocfd.fy, tocfd.fy + tocfd.fy_Len, it);
            std::advance(it, tocfd.fy_Len);
            std::copy(tocfd.fz, tocfd.fz + tocfd.fz_Len, it);
            std::advance(it, tocfd.fz_Len);

            std::copy(tocfd.pxForce, tocfd.pxForce + tocfd.pxForce_Len, it);
            std::advance(it, tocfd.pxForce_Len);
            std::copy(tocfd.pyForce, tocfd.pyForce + tocfd.pyForce_Len, it);
            std::advance(it, tocfd.pyForce_Len);
            std::copy(tocfd.pzForce, tocfd.pzForce + tocfd.pzForce_Len, it);
            std::advance(it, tocfd.pzForce_Len);

            // clang-format off
            std::copy(tocfd.pOrientation,
                      tocfd.pOrientation + tocfd.pOrientation_Len, it);
            // clang-format on
        }

        // Broadcast data to all influenced procs from the root process
        const auto& procs = data.info().procs;
        const int tag = 1001;
        if (data.info().is_root_proc) {
            BL_PROFILE(
                "amr-wind::actuator::ComputeForceOp<TurbineFast>::scatter2");
            for (const int ip : procs) {
                if (ip == data.info().root_proc) continue;

                amrex::ParallelDescriptor::Send(
                    buf.data(), dsize, ip, tag, data.meta().tcomm);
            }
        } else {
            BL_PROFILE(
                "amr-wind::actuator::ComputeForceOp<TurbineFast>::scatter2");
            amrex::ParallelDescriptor::Recv(
                buf.data(), dsize, data.info().root_proc, tag,
                data.meta().tcomm);
        }

        // Populate the actuator grid data structures with data from the MPI
        // send/recv buffer.
        {
            BL_PROFILE(
                "amr-wind::actuator::ComputeForceOp<TurbineFast>::scatter3");
            const auto& bp = data.info().base_pos;
            auto& grid = data.grid();
            const auto& npts = grid.pos.size();
            const auto& rho = data.meta().density;
            const size_t ifx = 0;
            const size_t ify = ifx + npts;
            const size_t ifz = ify + npts;
            const size_t ipx = ifz + npts;
            const size_t ipy = ipx + npts;
            const size_t ipz = ipy + npts;
            const size_t iori = ipz + npts;

            for (int i = 0; i < npts; ++i) {
                // Aerodynamic force vectors. Flip sign to get force on fluid.
                // Divide by density as the source term computation will
                // multiply by density before adding to momentum equation.
                //
                grid.force[i].x() =
                    -static_cast<amrex::Real>(buf[ifx + i]) / rho;
                grid.force[i].y() =
                    -static_cast<amrex::Real>(buf[ify + i]) / rho;
                grid.force[i].z() =
                    -static_cast<amrex::Real>(buf[ifz + i]) / rho;

                // Position vectors of the actuator nodes. Add shift to base
                // locations.
                grid.pos[i].x() =
                    static_cast<amrex::Real>(buf[ipx + i]) + bp.x();
                grid.pos[i].y() =
                    static_cast<amrex::Real>(buf[ipy + i]) + bp.y();
                grid.pos[i].z() =
                    static_cast<amrex::Real>(buf[ipz + i]) + bp.z();

                // Copy over the orientation matrix
                //
                // Note that we transpose the orientation matrix when copying
                // from OpenFAST to AMR-Wind Tensor data structure. This is done
                // so that post-multiplication of vector transforms from global
                // to local reference frame.
                const size_t off = iori + i * AMREX_SPACEDIM * AMREX_SPACEDIM;
                for (int j = 0; j < AMREX_SPACEDIM; ++j)
                    for (int k = 0; k < AMREX_SPACEDIM; ++k)
                        grid.orientation[i][j * AMREX_SPACEDIM + k] =
                            static_cast<amrex::Real>(
                                buf[off + j + k * AMREX_SPACEDIM]);
            }

            // Extract the rotor center of rotation
            auto& meta = data.meta();
            meta.rot_center = grid.pos[0];

            // Rotor non-rotating reference frame
            const auto xvec = grid.orientation[0].x().unit();
            const auto yvec = vs::Vector::khat() ^ xvec;
            const auto zvec = xvec ^ yvec;
            meta.rotor_frame.rows(xvec, yvec.unit(), zvec.unit());
        }
    }
};

template <typename SrcTrait>
struct ProcessOutputsOp<TurbineFast, SrcTrait>
{
private:
    typename TurbineFast::DataType& m_data;

    //! Path to the output directory (specified by Actuator physics class)
    std::string m_out_dir{""};

    //! NetCDF output filename for this turbine
    std::string m_nc_filename{""};

    //! Output frequency (specified in input file)
    int m_out_freq{10};

public:
    explicit ProcessOutputsOp(typename TurbineFast::DataType& data)
        : m_data(data)
    {}

    void read_io_options(const utils::ActParser& pp)
    {
        pp.query("output_frequency", m_out_freq);
    }

    void prepare_outputs(const std::string& out_dir)
    {
        m_nc_filename = out_dir + "/" + m_data.info().label + ".nc";
        utils::prepare_netcdf_file(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid());
    }

    void write_outputs()
    {
        const auto& time = m_data.sim().time();
        const int tidx = time.time_index();
        if ((m_out_freq > 0) && (tidx % m_out_freq != 0)) return;

        utils::write_netcdf(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid(),
            time.new_time());
    }
};

} // namespace ops
} // namespace actuator
} // namespace amr_wind

#endif /* TURBINE_FAST_OPS_H */
