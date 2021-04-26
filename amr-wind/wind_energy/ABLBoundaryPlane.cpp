#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABLBoundaryPlane.H"
#include "amr-wind/wind_energy/ABLFillInflow.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"

namespace amr_wind {

#ifdef AMR_WIND_USE_NETCDF
namespace {

//! Return closest index (from lower) of value in vector
AMREX_FORCE_INLINE int
closest_index(const amrex::Vector<amrex::Real>& vec, const amrex::Real value)
{
    auto const it = std::upper_bound(vec.begin(), vec.end(), value);
    AMREX_ALWAYS_ASSERT(it != vec.end());

    const int idx = std::distance(vec.begin(), it);
    return std::max(idx - 1, 0);
}

//! Return indices perpendicular to normal
template <typename T = amrex::GpuArray<int, 2>>
AMREX_FORCE_INLINE T perpendicular_idx(const int normal)
{
    switch (normal) {
    case 0:
        return T{1, 2};
    case 1:
        return T{0, 2};
    case 2:
        return T{0, 1};
    default:
        amrex::Abort("Invalid normal value to determine perpendicular indices");
    }
    return T{-1, -1};
}

//! Return offset vector
AMREX_FORCE_INLINE amrex::IntVect offset(const int face_dir, const int normal)
{
    amrex::IntVect offset(amrex::IntVect::TheDimensionVector(normal));
    if (face_dir == 1) {
        for (auto& o : offset) o *= -1;
    }
    return offset;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
plane_idx(const int i, const int j, const int k, const int perp, const int lo)
{
    return (static_cast<int>(perp == 0) * i + static_cast<int>(perp == 1) * j +
            static_cast<int>(perp == 2) * k) -
           lo;
}

AMREX_FORCE_INLINE std::string level_name(int lev)
{
    return "level_" + std::to_string(lev);
}
} // namespace

void InletData::resize(const int size)
{
    m_data_n.resize(size);
    m_data_np1.resize(size);
    m_data_interp.resize(size);
}

void InletData::define_plane(const amrex::Orientation ori)
{
    m_data_n[ori] = std::unique_ptr<PlaneVector>(new PlaneVector);
    m_data_np1[ori] = std::unique_ptr<PlaneVector>(new PlaneVector);
    m_data_interp[ori] = std::unique_ptr<PlaneVector>(new PlaneVector);
}

void InletData::define_level_data(
    const amrex::Orientation ori, const amrex::Box& bx, const size_t nc)
{
    if (not this->is_populated(ori)) return;
    m_data_n[ori]->push_back(amrex::FArrayBox(bx, nc));
    m_data_np1[ori]->push_back(amrex::FArrayBox(bx, nc));
    m_data_interp[ori]->push_back(amrex::FArrayBox(bx, nc));
}

void InletData::read_data(
    ncutils::NCGroup& grp,
    const amrex::Orientation ori,
    const int lev,
    const Field* fld,
    const amrex::Real time,
    const amrex::Vector<amrex::Real>& times)
{
    const size_t nc = fld->num_comp();
    const int nstart = m_components[fld->id()];

    const int idx = closest_index(times, time);
    const int idxp1 = idx + 1;
    m_tn = times[idx];
    m_tnp1 = times[idxp1];
    AMREX_ALWAYS_ASSERT(((m_tn <= time) and (time <= m_tnp1)));

    const int normal = ori.coordDir();
    const amrex::GpuArray<int, 2> perp = perpendicular_idx(normal);

    const auto& bx = (*m_data_n[ori])[lev].box();
    const auto& lo = bx.loVect();
    const size_t n0 = bx.length(perp[0]);
    const size_t n1 = bx.length(perp[1]);

    amrex::Vector<size_t> start{
        static_cast<size_t>(idx), static_cast<size_t>(lo[perp[0]]),
        static_cast<size_t>(lo[perp[1]]), 0};
    amrex::Vector<size_t> count{1, n0, n1, nc};
    amrex::Vector<amrex::Real> buffer(n0 * n1 * nc);
    grp.var(fld->name()).get(buffer.data(), start, count);

    const auto& datn = ((*m_data_n[ori])[lev]).array();
    auto d_buffer = buffer.dataPtr();
    amrex::LoopOnCpu(bx, nc, [=](int i, int j, int k, int n) noexcept {
        const int i0 = plane_idx(i, j, k, perp[0], lo[perp[0]]);
        const int i1 = plane_idx(i, j, k, perp[1], lo[perp[1]]);
        datn(i, j, k, n + nstart) = d_buffer[((i0 * n1) + i1) * nc + n];
    });

    start[0] = static_cast<size_t>(idxp1);
    grp.var(fld->name()).get(buffer.data(), start, count);

    const auto& datnp1 = ((*m_data_np1[ori])[lev]).array();
    amrex::LoopOnCpu(bx, nc, [=](int i, int j, int k, int n) noexcept {
        const int i0 = plane_idx(i, j, k, perp[0], lo[perp[0]]);
        const int i1 = plane_idx(i, j, k, perp[1], lo[perp[1]]);
        datnp1(i, j, k, n + nstart) = d_buffer[((i0 * n1) + i1) * nc + n];
    });

    ((*m_data_n[ori])[lev]).prefetchToDevice();
    ((*m_data_np1[ori])[lev]).prefetchToDevice();
}

void InletData::interpolate(const amrex::Real time)
{
    m_tinterp = time;
    for (amrex::OrientationIter oit; oit; ++oit) {
        auto ori = oit();
        if (not this->is_populated(ori)) continue;

        const int nlevels = m_data_n[ori]->size();
        for (int lev = 0; lev < nlevels; ++lev) {

            const auto& datn = (*m_data_n[ori])[lev];
            const auto& datnp1 = (*m_data_np1[ori])[lev];
            auto& dati = (*m_data_interp[ori])[lev];

            dati.linInterp<amrex::RunOn::Device>(
                datn, 0, datnp1, 0, m_tn, m_tnp1, m_tinterp, dati.box(), 0,
                dati.nComp());
        }
    }
}

bool InletData::is_populated(amrex::Orientation ori) const
{
    return m_data_n[ori] != nullptr;
}
#endif

ABLBoundaryPlane::ABLBoundaryPlane(CFDSim& sim)
    : m_time(sim.time()), m_repo(sim.repo()), m_mesh(sim.mesh())
{
    amrex::ParmParse pp("ABL");
    int pp_io_mode = -1;
    pp.query("bndry_io_mode", pp_io_mode);
    switch (pp_io_mode) {
    case 0:
        m_io_mode = io_mode::output;
        m_is_initialized = true;
        break;
    case 1:
        m_io_mode = io_mode::input;
        m_is_initialized = true;
        break;
    default:
        m_io_mode = io_mode::undefined;
        m_is_initialized = false;
        return;
    }

#ifndef AMR_WIND_USE_NETCDF
    if (m_is_initialized) {
        amrex::Abort(
            "ABLBoundaryPlane capability with IO mode requires NetCDF");
    }
#else
    pp.query("bndry_write_frequency", m_write_frequency);
    pp.queryarr("bndry_planes", m_planes);
    pp.query("bndry_output_start_time", m_out_start_time);
    pp.queryarr("bndry_var_names", m_var_names);
    pp.get("bndry_file", m_filename);

#endif
}

void ABLBoundaryPlane::post_init_actions()
{
    if (!m_is_initialized) return;
    initialize_data();
    write_header();
    write_file();
    read_header();
    read_file();
}

void ABLBoundaryPlane::pre_advance_work()
{
    if (!m_is_initialized) return;
    read_file();
}

void ABLBoundaryPlane::post_advance_work()
{
    if (!m_is_initialized) return;
    write_file();
}

void ABLBoundaryPlane::initialize_data()
{
#ifdef AMR_WIND_USE_NETCDF
    for (const auto& plane : m_planes) {
        amrex::Vector<std::string> valid_planes{"xlo", "ylo"};

        if ((std::find(valid_planes.begin(), valid_planes.end(), plane) ==
             valid_planes.end())) {
            throw std::runtime_error(
                "Requested plane (" + plane +
                ") does not exist. Pick one of [xlo, ylo].");
        }
    }

    for (const auto& fname : m_var_names) {
        if (m_repo.field_exists(fname)) {
            auto& fld = m_repo.get_field(fname);
            if (m_io_mode == io_mode::input) {
                fld.register_fill_patch_op<ABLFillInflow>(
                    m_mesh, m_time, *this);
            }
            m_fields.emplace_back(&fld);
        } else {
            amrex::Abort(
                "ABLBoundaryPlane: invalid variable requested: " + fname);
        }
    }
#endif
}

void ABLBoundaryPlane::write_header()
{
#ifdef AMR_WIND_USE_NETCDF
    if (m_io_mode != io_mode::output) return;

    amrex::Print() << "Creating output NetCDF file: " << m_filename
                   << std::endl;

    auto ncf = ncutils::NCFile::create_par(
        m_filename, NC_CLOBBER | NC_NETCDF4 | NC_MPIIO,
        amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

    ncf.enter_def_mode();
    ncf.def_dim("sdim", 1);
    ncf.def_dim("pdim", 2);
    ncf.def_dim("vdim", 3);
    ncf.def_dim("nt", NC_UNLIMITED);
    ncf.def_var("time", NC_DOUBLE, {"nt"});

    const int nlevels = m_repo.num_active_levels();
    for (amrex::OrientationIter oit; oit; ++oit) {
        auto ori = oit();
        const std::string plane = m_plane_names[ori];

        if (std::find(m_planes.begin(), m_planes.end(), plane) ==
            m_planes.end())
            continue;

        auto plane_grp = ncf.def_group(plane);

        const int normal = ori.coordDir();
        auto v_normal = plane_grp.def_scalar("normal", NC_INT);
        v_normal.put(&normal);

        const int face_dir = ori.faceDir();
        auto v_face = plane_grp.def_scalar("side", NC_INT);
        v_face.put(&face_dir);

        const amrex::Vector<int> perp =
            perpendicular_idx<amrex::Vector<int>>(normal);
        auto v_perp = plane_grp.def_var("perpendicular", NC_INT, {"pdim"});
        v_perp.put(perp.data());

        for (int lev = 0; lev < nlevels; ++lev) {

            auto lev_grp = plane_grp.def_group(level_name(lev));

            const amrex::Box& minBox = m_mesh.boxArray(lev).minimalBox();
            lev_grp.def_dim("nx", minBox.length(0));
            lev_grp.def_dim("ny", minBox.length(1));
            lev_grp.def_dim("nz", minBox.length(2));

            lev_grp.def_var("lengths", NC_DOUBLE, {"pdim"});
            lev_grp.def_var("lo", NC_DOUBLE, {"pdim"});
            lev_grp.def_var("hi", NC_DOUBLE, {"pdim"});
            lev_grp.def_var("dx", NC_DOUBLE, {"pdim"});

            const amrex::Vector<std::string> dirs{"nx", "ny", "nz"};
            for (auto* fld : m_fields) {
                const std::string name = fld->name();
                if (fld->num_comp() == 1) {
                    lev_grp.def_var(
                        name, NC_DOUBLE, {"nt", dirs[perp[0]], dirs[perp[1]]});
                } else if (fld->num_comp() == AMREX_SPACEDIM) {
                    lev_grp.def_var(
                        name, NC_DOUBLE,
                        {"nt", dirs[perp[0]], dirs[perp[1]], "vdim"});
                }
            }
        }
    }
    ncf.put_attr("title", m_title);
    ncf.exit_def_mode();

    // Populate coordinates
    for (auto& plane_grp : ncf.all_groups()) {
        int normal;
        plane_grp.var("normal").get(&normal);
        const amrex::GpuArray<int, 2> perp = perpendicular_idx(normal);

        for (int lev = 0; lev < nlevels; ++lev) {
            auto lev_grp = plane_grp.group(level_name(lev));

            const auto& dx = m_mesh.Geom(lev).CellSizeArray();
            const amrex::Box& minBox = m_mesh.boxArray(lev).minimalBox();
            const auto& lo = minBox.loVect();
            const auto& hi = minBox.hiVect();
            const amrex::Vector<amrex::Real> pdx{{dx[perp[0]], dx[perp[1]]}};
            const amrex::Vector<amrex::Real> los{
                {lo[perp[0]] * dx[perp[0]], lo[perp[1]] * dx[perp[1]]}};
            const amrex::Vector<amrex::Real> his{
                {(hi[perp[0]] + 1) * dx[perp[0]],
                 (hi[perp[1]] + 1) * dx[perp[1]]}};
            const amrex::Vector<amrex::Real> lengths{
                {minBox.length(perp[0]) * dx[perp[0]],
                 minBox.length(perp[1]) * dx[perp[1]]}};

            lev_grp.var("lengths").put(lengths.data());
            lev_grp.var("lo").put(los.data());
            lev_grp.var("hi").put(his.data());
            lev_grp.var("dx").put(pdx.data());
        }
    }

    amrex::Print() << "NetCDF file created successfully: " << m_filename
                   << std::endl;
#endif
}

void ABLBoundaryPlane::write_file()
{
#ifdef AMR_WIND_USE_NETCDF
    const amrex::Real time = m_time.new_time();
    const int t_step = m_time.time_index();

    // Only output data if at the desired timestep
    if ((t_step % m_write_frequency != 0) || ((m_io_mode != io_mode::output)) ||
        (time < m_out_start_time))
        return;

    amrex::Print() << "\nWriting NetCDF file " << m_filename << " at time "
                   << time << std::endl;

    auto ncf = ncutils::NCFile::open_par(
        m_filename, NC_WRITE | NC_NETCDF4 | NC_MPIIO,
        amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

    auto v_time = ncf.var("time");
    v_time.par_access(NC_COLLECTIVE);
    const size_t nt = ncf.dim("nt").len();
    v_time.put(&time, {nt}, {1});

    for (auto* fld : m_fields) {
        fld->fillpatch(m_time.current_time());
    }

    const int nlevels = m_repo.num_active_levels();
    for (amrex::OrientationIter oit; oit; ++oit) {
        auto ori = oit();
        const std::string plane = m_plane_names[ori];

        if (std::find(m_planes.begin(), m_planes.end(), plane) ==
            m_planes.end())
            continue;

        for (auto* fld : m_fields) {
            for (int lev = 0; lev < nlevels; ++lev) {
                auto grp = ncf.group(plane).group(level_name(lev));
                write_data(grp, ori, lev, fld);
            }
        }
    }

    m_out_counter++;
#endif
}

void ABLBoundaryPlane::read_header()
{
#ifdef AMR_WIND_USE_NETCDF
    if (m_io_mode != io_mode::input) return;

    // FIXME Do not support multi-level input mode yet.
    // this is due to interpolation issues at the coarse-fine interface
    if (m_repo.num_active_levels() > 1) {
        amrex::Abort("Not supporting multi-level input mode yet.");
    }

    amrex::Print() << "Reading input NetCDF file: " << m_filename << std::endl;
    auto ncf = ncutils::NCFile::open_par(
        m_filename, NC_NOWRITE | NC_NETCDF4 | NC_MPIIO,
        amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

    // Store the input file times and reset to start at 0
    const size_t nt = ncf.dim("nt").len();
    m_in_times.resize(nt);
    ncf.var("time").get(m_in_times.data());

    // Sanity check the input file time
    AMREX_ALWAYS_ASSERT(m_in_times[0] <= m_time.current_time());

    const int nlevels = m_repo.num_active_levels();
    m_in_data.resize(6);
    for (auto& plane_grp : ncf.all_groups()) {
        int normal, face_dir;
        plane_grp.var("normal").get(&normal);
        plane_grp.var("side").get(&face_dir);
        const amrex::GpuArray<int, 2> perp = perpendicular_idx(normal);
        const amrex::Orientation ori(
            normal, amrex::Orientation::Side(face_dir));

        m_in_data.define_plane(ori);
        for (int lev = 0; lev < nlevels; ++lev) {
            auto lev_grp = plane_grp.group(level_name(lev));

            // sanity checks to ensure grid-to-grid matching
            const amrex::Box& minBox = m_mesh.boxArray(lev).minimalBox();
            const auto& lo = minBox.loVect();
            const auto& hi = minBox.hiVect();
            const auto& dx = m_mesh.Geom(lev).CellSizeArray();
            const amrex::Vector<amrex::Real> pdx{{dx[perp[0]], dx[perp[1]]}};
            const amrex::Vector<amrex::Real> los{
                {lo[perp[0]] * pdx[0], lo[perp[1]] * pdx[1]}};
            const amrex::Vector<amrex::Real> his{
                {(hi[perp[0]] + 1) * pdx[0], (hi[perp[1]] + 1) * pdx[1]}};
            const amrex::Vector<amrex::Real> lengths{
                {minBox.length(perp[0]) * pdx[0],
                 minBox.length(perp[1]) * pdx[1]}};

            amrex::Vector<amrex::Real> nc_dat{{0, 0}};
            lev_grp.var("lengths").get(nc_dat.data());
            AMREX_ALWAYS_ASSERT(nc_dat == lengths);
            lev_grp.var("lo").get(nc_dat.data());
            AMREX_ALWAYS_ASSERT(nc_dat == los);
            lev_grp.var("hi").get(nc_dat.data());
            AMREX_ALWAYS_ASSERT(nc_dat == his);
            lev_grp.var("dx").get(nc_dat.data());
            AMREX_ALWAYS_ASSERT(nc_dat == pdx);

            // Create the data structures for the input data
            amrex::IntVect plo(lo);
            amrex::IntVect phi(hi);
            plo[normal] = ori.isHigh() ? lo[normal] + 1 : -1;
            phi[normal] = ori.isHigh() ? hi[normal] + 1 : -1;
            const amrex::Box pbx(plo, phi);
            size_t nc = 0;
            for (auto* fld : m_fields) {
                m_in_data.component(fld->id()) = nc;
                nc += fld->num_comp();
            }
            m_in_data.define_level_data(ori, pbx, nc);
        }
    }

    amrex::Print() << "NetCDF file read successfully: " << m_filename
                   << std::endl;
#endif
}

void ABLBoundaryPlane::read_file()
{
#ifdef AMR_WIND_USE_NETCDF
    if (m_io_mode != io_mode::input) return;

    // populate planes and interpolate
    const amrex::Real time = m_time.new_time();
    AMREX_ALWAYS_ASSERT((m_in_times[0] <= time) and (time < m_in_times.back()));

    if (!((m_in_data.tn() <= time) && (time < m_in_data.tnp1()))) {

        auto ncf = ncutils::NCFile::open_par(
            m_filename, NC_NOWRITE | NC_NETCDF4 | NC_MPIIO,
            amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

        const int nlevels = m_repo.num_active_levels();
        for (amrex::OrientationIter oit; oit; ++oit) {
            auto ori = oit();
            if (not m_in_data.is_populated(ori)) continue;

            const std::string plane = m_plane_names[ori];
            for (auto* fld : m_fields) {
                for (int lev = 0; lev < nlevels; ++lev) {
                    auto grp = ncf.group(plane).group(level_name(lev));
                    m_in_data.read_data(grp, ori, lev, fld, time, m_in_times);
                }
            }
        }
    }

    m_in_data.interpolate(time);
#endif
}

void ABLBoundaryPlane::populate_data(
    const int lev,
    const amrex::Real time,
    Field& fld,
    amrex::MultiFab& mfab) const
{
#ifdef AMR_WIND_USE_NETCDF

    if (m_io_mode != io_mode::input) return;

    AMREX_ALWAYS_ASSERT(
        ((m_in_data.tn() <= time) || (time <= m_in_data.tnp1())));
    AMREX_ALWAYS_ASSERT(amrex::Math::abs(time - m_in_data.tinterp()) < 1e-12);

    for (amrex::OrientationIter oit; oit; ++oit) {
        auto ori = oit();
        if ((not m_in_data.is_populated(ori)) or
            (fld.bc_type()[ori] != BC::mass_inflow))
            continue;

        const int normal = ori.coordDir();
        const amrex::GpuArray<int, 2> perp = perpendicular_idx(normal);

        const size_t nc = mfab.nComp();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(mfab, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {

            const auto& sbx = mfi.growntilebox(1);
            auto& src = m_in_data.interpolate_data(ori, lev);
            const auto& bx = sbx & src.box();
            if (bx.isEmpty()) continue;

            const auto& dest = mfab.array(mfi);
            const auto& src_arr = src.array();
            const int nstart = m_in_data.component(fld.id());
            amrex::ParallelFor(
                bx, nc,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    dest(i, j, k, n) = src_arr(i, j, k, n + nstart);
                });

            // Fill the edges
            const auto& lo = bx.loVect();
            const auto& hi = bx.hiVect();

            // For xlo/ylo combination (currently the only valid
            // combination), this is perp[0] (FIXME for future)
            const int pp = perp[0];

            {
                amrex::IntVect elo(lo);
                amrex::IntVect ehi(hi);
                ehi[pp] = lo[pp];
                const amrex::Box ebx(elo, ehi);

                amrex::GpuArray<int, 3> v_offset{{pp == 0, pp == 1, pp == 2}};
                amrex::ParallelFor(
                    ebx, nc,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                        dest(
                            i - v_offset[0], j - v_offset[1], k - v_offset[2],
                            n) = dest(i, j, k, n);
                    });
            }

            {
                amrex::IntVect elo(lo);
                amrex::IntVect ehi(hi);
                elo[pp] = hi[pp];
                const amrex::Box ebx(elo, ehi);

                amrex::GpuArray<int, 3> v_offset{{pp == 0, pp == 1, pp == 2}};
                amrex::ParallelFor(
                    ebx, nc,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                        dest(
                            i + v_offset[0], j + v_offset[1], k + v_offset[2],
                            n) = dest(i, j, k, n);
                    });
            }
        }
    }
#else
    amrex::ignore_unused(lev, time, fld, mfab);
#endif
}

#ifdef AMR_WIND_USE_NETCDF
void ABLBoundaryPlane::write_data(
    ncutils::NCGroup& grp,
    const amrex::Orientation ori,
    const int lev,
    const Field* fld)
{
    // Plane info
    const int normal = ori.coordDir();
    const amrex::GpuArray<int, 2> perp = perpendicular_idx(normal);
    const amrex::IntVect v_offset = offset(ori.faceDir(), normal);

    // Field info
    const std::string name = fld->name();
    const size_t nc = fld->num_comp();

    // Domain info
    const amrex::Box& domain = m_mesh.Geom(lev).Domain();
    const auto& dlo = domain.loVect();
    AMREX_ALWAYS_ASSERT(dlo[0] == 0 and dlo[1] == 0 and dlo[2] == 0);

    grp.var(name).par_access(NC_COLLECTIVE);

    // FIXME optimization
    // - move buffer outside this function, probably best as a member
    // - place in object to access as ori/lev/fld
    // - sizing and start/counts should be done only on init and regrid
    const int n_buffers = m_mesh.boxArray(lev).size();
    amrex::Vector<BufferData> buffers(n_buffers);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi((*fld)(lev), false); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.tilebox();
        const auto& blo = bx.loVect();
        const auto& bhi = bx.hiVect();

        if (blo[normal] == dlo[normal]) {
            amrex::IntVect lo(blo);
            amrex::IntVect hi(bhi);
            lo[normal] = dlo[normal];
            hi[normal] = dlo[normal];
            const amrex::Box lbx(lo, hi);

            const size_t n0 = hi[perp[0]] - lo[perp[0]] + 1;
            const size_t n1 = hi[perp[1]] - lo[perp[1]] + 1;

            auto& buffer = buffers[mfi.index()];
            buffer.data.resize(n0 * n1 * nc);

            auto const& fld_arr = (*fld)(lev).array(mfi);
            impl_buffer_field(
                lbx, n1, nc, perp, v_offset, fld_arr, buffer.data);
            amrex::Gpu::streamSynchronize();

            buffer.start = {
                m_out_counter, static_cast<size_t>(lo[perp[0]]),
                static_cast<size_t>(lo[perp[1]]), 0};
            buffer.count = {1, n0, n1, nc};
        }
    }

    for (const auto& buffer : buffers) {
        grp.var(name).put(buffer.data.dataPtr(), buffer.start, buffer.count);
    }
}

void ABLBoundaryPlane::impl_buffer_field(
    const amrex::Box& bx,
    const int n1,
    const int nc,
    const amrex::GpuArray<int, 2>& perp,
    const amrex::IntVect& v_offset,
    const amrex::Array4<const amrex::Real>& fld,
    amrex::Gpu::ManagedVector<amrex::Real>& buffer)
{
    auto d_buffer = buffer.dataPtr();
    const auto lo = bx.loVect3d();
    amrex::ParallelFor(
        bx, nc, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            const int i0 = plane_idx(i, j, k, perp[0], lo[perp[0]]);
            const int i1 = plane_idx(i, j, k, perp[1], lo[perp[1]]);
            d_buffer[((i0 * n1) + i1) * nc + n] =
                0.5 * (fld(i, j, k, n) + fld(i - v_offset[0], j - v_offset[1],
                                             k - v_offset[2], n));
        });
}
#endif
} // namespace amr_wind
