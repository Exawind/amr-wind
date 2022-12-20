#ifndef HOS_WAVES_OPS_H
#define HOS_WAVES_OPS_H

#include "amr-wind/ocean_waves/relaxation_zones/HOSWaves.H"
#include "amr-wind/ocean_waves/OceanWavesOps.H"
#include "amr-wind/ocean_waves/OceanWaves.H"
#include "amr-wind/ocean_waves/relaxation_zones/relaxation_zones_ops.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/equation_systems/BCOps.H"

namespace amr_wind::ocean_waves::ops {

void ReadHOSFileLev(
    int lev,
    const std::string& HOS_prefix,
    int HOS_n,
    amrex::Real& HOS_t,
    amrex::Real& HOS_dt,
    RealList& eta,
    RealList& u,
    RealList& v,
    RealList& w,
    int& HOS_nx,
    amrex::Real& HOS_Lx,
    int& HOS_ny,
    amrex::Real& HOS_Ly,
    int& HOS_nz,
    amrex::Real& HOS_zmin,
    amrex::Real& HOS_zmax)
{
    // Convert file string to stream
    std::stringstream fname;
    fname << HOS_prefix << "_lev" << lev << "_" << HOS_n << ".txt";
    std::ifstream is(fname.str());
    // Read metadata from file
    std::string tmp;
    // Get initial time
    std::getline(is, tmp, '=');
    std::getline(is, tmp);
    HOS_t = std::stof(tmp);
    // Get dt
    std::getline(is, tmp, '=');
    std::getline(is, tmp);
    HOS_dt = std::stof(tmp);
    // Get nx, Lx
    std::getline(is, tmp, '=');
    std::getline(is, tmp, ',');
    HOS_nx = std::stoi(tmp);
    std::getline(is, tmp, '=');
    std::getline(is, tmp);
    HOS_Lx = std::stof(tmp);
    // Get ny, Ly
    std::getline(is, tmp, '=');
    std::getline(is, tmp, ',');
    HOS_ny = std::stoi(tmp);
    std::getline(is, tmp, '=');
    std::getline(is, tmp);
    HOS_Ly = std::stof(tmp);
    // Get nz, zmin, zmax
    std::getline(is, tmp, '=');
    std::getline(is, tmp, ',');
    HOS_nz = std::stoi(tmp);
    std::getline(is, tmp, '=');
    std::getline(is, tmp, ',');
    HOS_zmin = std::stof(tmp);
    std::getline(is, tmp, '=');
    std::getline(is, tmp);
    HOS_zmax = std::stof(tmp);

    // Modify Lx, Ly to account for node - cell difference
    /*HOS_Lx *= (((amrex::Real)HOS_nx) / (amrex::Real)(HOS_nx -
    std::pow(2,lev)));
    HOS_Ly *= (((amrex::Real)HOS_ny) / (amrex::Real)(HOS_ny -
    std::pow(2,lev)));*/

    // Allocate arrays for storage
    eta.resize(static_cast<long>(HOS_nx) * HOS_ny);
    u.resize(static_cast<long>(HOS_nx) * HOS_ny * HOS_nz);
    v.resize(static_cast<long>(HOS_nx) * HOS_ny * HOS_nz);
    w.resize(static_cast<long>(HOS_nx) * HOS_ny * HOS_nz);
    // Skip key
    std::getline(is, tmp);
    // Read interface heights and velocities
    for (int ilat = 0; ilat < HOS_nx * HOS_ny; ++ilat) {
        // Get eta for current point
        is >> eta[ilat];
        // Get u, v, w for full depth of 2D point
        for (int ivert = 0; ivert < HOS_nz; ++ivert) {
            is >> u[ilat * HOS_nz + ivert] >> v[ilat * HOS_nz + ivert] >>
                w[ilat * HOS_nz + ivert];
        }
    }
}

void ReadHOSFileLev(
    int lev,
    const std::string& HOS_prefix,
    int HOS_n,
    RealList& eta,
    RealList& u,
    RealList& v,
    RealList& w,
    int& HOS_nx,
    amrex::Real& HOS_Lx,
    int& HOS_ny,
    amrex::Real& HOS_Ly,
    int& HOS_nz,
    amrex::Real& HOS_zmin,
    amrex::Real& HOS_zmax)
{
    amrex::Real tmp1, tmp2;
    ReadHOSFileLev(
        lev, HOS_prefix, HOS_n, tmp1, tmp2, eta, u, v, w, HOS_nx, HOS_Lx,
        HOS_ny, HOS_Ly, HOS_nz, HOS_zmin, HOS_zmax);
}

void StoreHOSDataLoop(
    HOSWaves::MetaType& wdata,
    amrex::Array4<amrex::Real> const& phi,
    amrex::Array4<amrex::Real> const& vel,
    const amrex::Real* dev_eta_data,
    const amrex::Real* dev_u_data,
    const amrex::Real* dev_v_data,
    const amrex::Real* dev_w_data,
    const amrex::Real HOS_Lx,
    const int HOS_nx,
    const amrex::Real HOS_Ly,
    const int HOS_ny,
    const amrex::Real HOS_zmin,
    const amrex::Real HOS_zmax,
    const int HOS_nz,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> problo,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
    const amrex::Box& vbx)
{
    const amrex::Real zsl = wdata.zsl;
    const amrex::Real HOS_dx = HOS_Lx / HOS_nx;
    const amrex::Real HOS_dy = HOS_Ly / HOS_ny;
    const amrex::Real HOS_dz = (HOS_zmax - HOS_zmin) / HOS_nz;
    const auto& gbx = grow(vbx, 3);

    amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real x = problo[0] + (i + 0.5) * dx[0];
        amrex::Real y = problo[1] + (j + 0.5) * dx[1];
        amrex::Real z = problo[2] + (k + 0.5) * dx[2];
        // Adjust x,y for periodicity of HOS domain
        if (x > HOS_Lx) {
            x -= HOS_Lx;
        }
        if (x < 0.0) {
            x += HOS_Lx;
        }
        if (y > HOS_Ly) {
            y -= HOS_Ly;
        }
        if (y < 0.0) {
            y += HOS_Ly;
        }
        // Determine lateral indices for interpolation
        // It's best to have the grids line up directly
        int ii1 = 0;
        amrex::Real HOS_x = HOS_dx * (ii1 + 0.5);
        while (HOS_x < x) {
            ++ii1;
            HOS_x = HOS_dx * (ii1 + 0.5);
        }
        int jj1 = 0;
        amrex::Real HOS_y = HOS_dy * (jj1 + 0.5);
        while (HOS_y < y) {
            ++jj1;
            HOS_y = HOS_dy * (jj1 + 0.5);
        }

        // Do interpolation for x, y and get eta
        // Shift indices left (because of the finding function)
        int ii = ii1 - 1;
        int jj = jj1 - 1;
        // Reset locations
        HOS_x = HOS_dx * (ii + 0.5);
        HOS_y = HOS_dy * (jj + 0.5);
        // Modify indices to avoid going out of bounds,
        // according to periodicity (assumed about HOS data)

        if (ii < 0) {
            ii += HOS_nx;
        }
        if (ii1 < 0) {
            ii1 += HOS_nx;
        }
        if (jj < 0) {
            jj += HOS_ny;
        }
        if (jj1 < 0) {
            jj1 += HOS_ny;
        }
        if (ii >= HOS_nx) {
            ii -= HOS_nx;
        }
        if (ii1 >= HOS_nx) {
            ii1 -= HOS_nx;
        }
        if (jj >= HOS_ny) {
            jj -= HOS_ny;
        }
        if (jj1 >= HOS_ny) {
            jj1 -= HOS_ny;
        }

        // Interpolation factors
        amrex::Real wx_hi = (x - HOS_x) / HOS_dx;
        amrex::Real wy_hi = (y - HOS_y) / HOS_dy;
        amrex::Real wx_lo = 1.0 - wx_hi;
        amrex::Real wy_lo = 1.0 - wy_hi;

        // Get data
        const amrex::Real eta_00 = dev_eta_data[jj + ii * HOS_ny];
        const amrex::Real eta_01 = dev_eta_data[jj1 + ii * HOS_ny];
        const amrex::Real eta_10 = dev_eta_data[jj + ii1 * HOS_ny];
        const amrex::Real eta_11 = dev_eta_data[jj1 + ii1 * HOS_ny];
        const amrex::Real eta_ =
            wx_lo * wy_lo * eta_00 + wx_lo * wy_hi * eta_01 +
            wx_hi * wy_lo * eta_10 + wx_hi * wy_hi * eta_11;

        phi(i, j, k) = eta_ + zsl - z;

        // Determine z location for interpolation based on phase
        if (phi(i, j, k) + 0.5 * dx[2] < 0) {
            // Fully gas (exceed range so that velocity = 0)
            z = HOS_zmax + 1;
        } else if (amrex::Math::abs(phi(i, j, k)) - 0.5 * dx[2] < 0) {
            // Mixed phase (assume liquid is on bottom)
            z -= 0.5 * phi(i, j, k);
        } // else: Fully liquid, z is unchanged

        if (z > HOS_zmax || z < HOS_zmin) {
            // Assume velocities are zero outside range
            vel(i, j, k, 0) = 0.;
            vel(i, j, k, 1) = 0.;
            vel(i, j, k, 2) = 0.;
        } else {
            // Determine vertical indices
            int kk = 0;
            amrex::Real HOS_z = HOS_zmin + HOS_dz * (kk + 0.5);
            while (HOS_z < z) {
                ++kk;
                HOS_z = HOS_zmin + HOS_dz * (kk + 0.5);
            }
            // Shift indices left
            kk = kk - 1;
            // Modify indices to avoid going out of bounds
            kk = amrex::max(0, amrex::min(HOS_nz - 2, kk));
            // Reset location
            HOS_z = HOS_zmin + HOS_dz * (kk + 0.5);
            // Interpolation factors
            amrex::Real wz_hi = (z - HOS_z) / HOS_dz;
            amrex::Real wz_lo = 1.0 - wz_hi;

            // Set up indices for clarity
            const int i000 = (jj + ii * HOS_ny) * HOS_nz + kk;
            const int i001 = (jj + ii * HOS_ny) * HOS_nz + kk + 1;
            const int i010 = (jj1 + ii * HOS_ny) * HOS_nz + kk;
            const int i100 = (jj + ii1 * HOS_ny) * HOS_nz + kk;
            const int i011 = (jj1 + ii * HOS_ny) * HOS_nz + kk + 1;
            const int i101 = (jj + ii1 * HOS_ny) * HOS_nz + kk + 1;
            const int i110 = (jj1 + ii1 * HOS_ny) * HOS_nz + kk;
            const int i111 = (jj1 + ii1 * HOS_ny) * HOS_nz + kk + 1;
            vel(i, j, k, 0) = wx_lo * wy_lo * wz_lo * dev_u_data[i000] +
                              wx_lo * wy_lo * wz_hi * dev_u_data[i001] +
                              wx_lo * wy_hi * wz_lo * dev_u_data[i010] +
                              wx_lo * wy_hi * wz_hi * dev_u_data[i011] +
                              wx_hi * wy_lo * wz_lo * dev_u_data[i100] +
                              wx_hi * wy_lo * wz_hi * dev_u_data[i101] +
                              wx_hi * wy_hi * wz_lo * dev_u_data[i110] +
                              wx_hi * wy_hi * wz_hi * dev_u_data[i111];
            vel(i, j, k, 1) = wx_lo * wy_lo * wz_lo * dev_v_data[i000] +
                              wx_lo * wy_lo * wz_hi * dev_v_data[i001] +
                              wx_lo * wy_hi * wz_lo * dev_v_data[i010] +
                              wx_lo * wy_hi * wz_hi * dev_v_data[i011] +
                              wx_hi * wy_lo * wz_lo * dev_v_data[i100] +
                              wx_hi * wy_lo * wz_hi * dev_v_data[i101] +
                              wx_hi * wy_hi * wz_lo * dev_v_data[i110] +
                              wx_hi * wy_hi * wz_hi * dev_v_data[i111];
            vel(i, j, k, 2) = wx_lo * wy_lo * wz_lo * dev_w_data[i000] +
                              wx_lo * wy_lo * wz_hi * dev_w_data[i001] +
                              wx_lo * wy_hi * wz_lo * dev_w_data[i010] +
                              wx_lo * wy_hi * wz_hi * dev_w_data[i011] +
                              wx_hi * wy_lo * wz_lo * dev_w_data[i100] +
                              wx_hi * wy_lo * wz_hi * dev_w_data[i101] +
                              wx_hi * wy_hi * wz_lo * dev_w_data[i110] +
                              wx_hi * wy_hi * wz_hi * dev_w_data[i111];

            // Set velocity to zero if no liquid present or if wave
            // initialization not desired
            if (phi(i, j, k) + 0.5 * dx[2] < 0.0) {
                vel(i, j, k, 0) = 0.0;
                vel(i, j, k, 1) = 0.0;
                vel(i, j, k, 2) = 0.0;
            }
        }
    });
}

template <>
struct ReadInputsOp<HOSWaves>
{
    void operator()(
        HOSWaves::DataType& data, const ::amr_wind::utils::MultiParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();
        relaxation_zones::read_inputs(wdata, info, pp);

        pp.get("HOS_files_prefix", wdata.HOS_prefix);
        pp.query("HOS_init_timestep", wdata.HOS_n0);
        wdata.HOS_n = wdata.HOS_n0;

        // Declare fields for HOS
        auto& hos_levelset =
            data.sim().repo().declare_field("hos_levelset", 1, 3, 1);
        auto& hos_velocity = data.sim().repo().declare_field(
            "hos_velocity", AMREX_SPACEDIM, 3, 1);

        // Extrapolation can work well when finer data is available
        hos_levelset.set_default_fillpatch_bc(data.sim().time());
        hos_velocity.set_default_fillpatch_bc(data.sim().time());
        // No operation works better when only coarse data is used
        /*hos_levelset.register_fill_patch_op<FieldFillPatchOps<FieldBCNoOp>>(
            data.sim().repo().mesh(), data.sim().time(), 0);
        hos_velocity.register_fill_patch_op<FieldFillPatchOps<FieldBCNoOp>>(
            data.sim().repo().mesh(), data.sim().time(), 0);*/
    }
}; // namespace ops

template <>
struct InitDataOp<HOSWaves>
{
    void
    operator()(HOSWaves::DataType& data, int level, const amrex::Geometry& geom)
    {
        auto& wdata = data.meta();
        bool init_waves = wdata.init_wave_field;

        auto& sim = data.sim();
        // cppcheck-suppress constVariable

        auto& m_levelset = sim.repo().get_field("levelset");
        auto& m_velocity = sim.repo().get_field("velocity");
        // Set up variables that are re-written at each level
        RealList eta, u, v, w;
        int HOS_nx, HOS_ny, HOS_nz;
        amrex::Real HOS_Lx, HOS_Ly, HOS_zmin, HOS_zmax;

        const auto& problo = geom.ProbLoArray();
        const auto& probhi = geom.ProbHiArray();
        const auto& dx = geom.CellSizeArray();
        // Read HOS data at current level
        ReadHOSFileLev(
            level, wdata.HOS_prefix, wdata.HOS_n, wdata.HOS_t, wdata.HOS_dt,
            eta, u, v, w, HOS_nx, HOS_Lx, HOS_ny, HOS_Ly, HOS_nz, HOS_zmin,
            HOS_zmax);

        // Check if current dimensions are compatible
        if (problo[0] < -1e-6 || probhi[0] > HOS_Lx * (1.0 + 1e-6) ||
            problo[1] < -1e-6 || probhi[1] > HOS_Ly * (1.0 + 1e-6)) {
            amrex::Abort(
                "HOS OceanWaves: Lateral dimensions incompatible, level " +
                std::to_string(level));
        }
        amrex::Gpu::DeviceVector<amrex::Real> dev_eta, dev_u, dev_v, dev_w;
        dev_eta.resize(static_cast<long>(HOS_nx) * HOS_ny);
        dev_u.resize(static_cast<long>(HOS_nx) * HOS_ny * HOS_nz);
        dev_v.resize(static_cast<long>(HOS_nx) * HOS_ny * HOS_nz);
        dev_w.resize(static_cast<long>(HOS_nx) * HOS_ny * HOS_nz);

        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, eta.begin(), eta.end(), dev_eta.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, u.begin(), u.end(), dev_u.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, v.begin(), v.end(), dev_v.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, w.begin(), w.end(), dev_w.begin());

        const amrex::Real* dev_eta_data = dev_eta.data();
        const amrex::Real* dev_u_data = dev_u.data();
        const amrex::Real* dev_v_data = dev_v.data();
        const amrex::Real* dev_w_data = dev_w.data();
        // Loop to populate field data
        for (amrex::MFIter mfi(m_levelset(level)); mfi.isValid(); ++mfi) {
            auto phi = m_levelset(level).array(mfi);
            auto vel = m_velocity(level).array(mfi);
            const auto& vbx = mfi.validbox();
            if (init_waves) {
                StoreHOSDataLoop(
                    wdata, phi, vel, dev_eta_data, dev_u_data, dev_v_data,
                    dev_w_data, HOS_Lx, HOS_nx, HOS_Ly, HOS_ny, HOS_zmin,
                    HOS_zmax, HOS_nz, problo, dx, vbx);
            } else {
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
        }
        // Average down from top level for consistency
        const int nlevels = sim.repo().num_active_levels();
        if (level == nlevels - 1) {
            for (int lev = nlevels - 1; lev > 0; --lev) {
                amrex::average_down(
                    m_velocity(lev), m_velocity(lev - 1), 0, AMREX_SPACEDIM,
                    sim.mesh().refRatio(lev - 1));
                amrex::average_down(
                    m_levelset(lev), m_levelset(lev - 1), 0, 1,
                    sim.mesh().refRatio(lev - 1));
            }
        }
    }
}; // namespace ocean_waves

template <>
struct UpdateRelaxZonesOp<HOSWaves>
{
    void operator()(HOSWaves::DataType& data)
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
        auto& hos_levelset = sim.repo().get_field("hos_levelset");
        auto& hos_velocity = sim.repo().get_field("hos_velocity");

        // Get value for time interpolation
        amrex::Real t_last = wdata.t_last;

        // Check if HOS time needs to be incremented
        bool read_flag = false;
        while (time > wdata.HOS_t) {
            // Increment time step
            ++wdata.HOS_n;
            // Increment time
            wdata.HOS_t += wdata.HOS_dt;
            // Flag to store new data
            read_flag = true;
        }
        // Signifies initialization from scratch or a restart
        if (wdata.t_last < 0.0) {
            read_flag = true;
            // Need to read in previous or current time for ow state
            while (time < wdata.HOS_t) {
                // Go back
                --wdata.HOS_n;
                wdata.HOS_t -= wdata.HOS_dt;
            }

            // Store the correct time for next interpolation
            wdata.t_last = wdata.HOS_t;
            // Modify time to get initial ow values correct
            time = wdata.HOS_t;

            // Initialize ow fields to 0 to be safe for time interp
            m_ow_levelset.setVal(0.0);
            m_ow_velocity.setVal(0.0);
        } else {
            // Record latest time as 'last' for next timestep
            wdata.t_last = time;
        }

        auto nlevels = sim.repo().num_active_levels();
        auto geom = sim.mesh().Geom();

        // Read HOS data if necessary
        if (read_flag) {
            // Set up variables that are re-written at each level
            RealList eta, u, v, w;
            int HOS_nx, HOS_ny, HOS_nz;
            amrex::Real HOS_Lx, HOS_Ly, HOS_zmin, HOS_zmax;
            for (int lev = 0; lev < nlevels; ++lev) {
                const auto& problo = geom[lev].ProbLoArray();
                const auto& dx = geom[lev].CellSizeArray();
                // Read HOS data at current level
                ReadHOSFileLev(
                    lev, wdata.HOS_prefix, wdata.HOS_n, eta, u, v, w, HOS_nx,
                    HOS_Lx, HOS_ny, HOS_Ly, HOS_nz, HOS_zmin, HOS_zmax);
                amrex::Gpu::DeviceVector<amrex::Real> dev_eta, dev_u, dev_v,
                    dev_w;
                dev_eta.resize(static_cast<long>(HOS_nx) * HOS_ny);
                dev_u.resize(static_cast<long>(HOS_nx) * HOS_ny * HOS_nz);
                dev_v.resize(static_cast<long>(HOS_nx) * HOS_ny * HOS_nz);
                dev_w.resize(static_cast<long>(HOS_nx) * HOS_ny * HOS_nz);

                amrex::Gpu::copy(
                    amrex::Gpu::hostToDevice, eta.begin(), eta.end(),
                    dev_eta.begin());
                amrex::Gpu::copy(
                    amrex::Gpu::hostToDevice, u.begin(), u.end(),
                    dev_u.begin());
                amrex::Gpu::copy(
                    amrex::Gpu::hostToDevice, v.begin(), v.end(),
                    dev_v.begin());
                amrex::Gpu::copy(
                    amrex::Gpu::hostToDevice, w.begin(), w.end(),
                    dev_w.begin());
                // Loop through multifab to interpolate and store in HOS fields
                for (amrex::MFIter mfi(m_ow_levelset(lev)); mfi.isValid();
                     ++mfi) {
                    auto HOS_phi = hos_levelset(lev).array(mfi);
                    auto HOS_vel = hos_velocity(lev).array(mfi);
                    const auto& vbx = mfi.validbox();

                    const amrex::Real* dev_eta_data = dev_eta.data();
                    const amrex::Real* dev_u_data = dev_u.data();
                    const amrex::Real* dev_v_data = dev_v.data();
                    const amrex::Real* dev_w_data = dev_w.data();

                    StoreHOSDataLoop(
                        wdata, HOS_phi, HOS_vel, dev_eta_data, dev_u_data,
                        dev_v_data, dev_w_data, HOS_Lx, HOS_nx, HOS_Ly, HOS_ny,
                        HOS_zmin, HOS_zmax, HOS_nz, problo, dx, vbx);
                }
            }
            // Average down to get fine information on coarse grid where
            // possible
            for (int lev = nlevels - 1; lev > 0; --lev) {
                amrex::average_down(
                    hos_velocity(lev), hos_velocity(lev - 1), 0, AMREX_SPACEDIM,
                    sim.mesh().refRatio(lev - 1));
                amrex::average_down(
                    hos_levelset(lev), hos_levelset(lev - 1), 0, 1,
                    sim.mesh().refRatio(lev - 1));
            }
            // Fill in across internal boundaries
            hos_velocity.fillpatch(0.0);
            hos_levelset.fillpatch(0.0);
        }

        // Temporally interpolate at every timestep to get target solution
        for (int lev = 0; lev < nlevels; ++lev) {
            for (amrex::MFIter mfi(m_ow_levelset(lev)); mfi.isValid(); ++mfi) {
                auto phi = m_ow_levelset(lev).array(mfi);
                auto vel = m_ow_velocity(lev).array(mfi);
                auto HOS_phi = hos_levelset(lev).array(mfi);
                auto HOS_vel = hos_velocity(lev).array(mfi);

                const amrex::Real HOS_t = wdata.HOS_t;
                const auto& gbx = mfi.growntilebox(3);
                amrex::ParallelFor(
                    gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        // Interpolate temporally every time
                        phi(i, j, k) += (HOS_phi(i, j, k) - phi(i, j, k)) *
                                        (time - t_last) /
                                        (HOS_t - t_last + 1e-16);
                        vel(i, j, k, 0) +=
                            (HOS_vel(i, j, k, 0) - vel(i, j, k, 0)) *
                            (time - t_last) / (HOS_t - t_last + 1e-16);
                        vel(i, j, k, 1) +=
                            (HOS_vel(i, j, k, 1) - vel(i, j, k, 1)) *
                            (time - t_last) / (HOS_t - t_last + 1e-16);
                        vel(i, j, k, 2) +=
                            (HOS_vel(i, j, k, 2) - vel(i, j, k, 2)) *
                            (time - t_last) / (HOS_t - t_last + 1e-16);
                    });
            }
        }
    }
};

} // namespace amr_wind::ocean_waves::ops

#endif /* HOS_WAVES_OPS_H */
