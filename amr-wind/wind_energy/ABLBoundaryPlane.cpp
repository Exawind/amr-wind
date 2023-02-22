#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABLBoundaryPlane.H"
#include "amr-wind/wind_energy/ABLFillInflow.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include <AMReX_PlotFileUtil.H>

namespace amr_wind {

namespace {

//! Return closest index (from lower) of value in vector
AMREX_FORCE_INLINE int
closest_index(const amrex::Vector<amrex::Real>& vec, const amrex::Real value)
{
    auto const it = std::upper_bound(vec.begin(), vec.end(), value);
    AMREX_ALWAYS_ASSERT(it != vec.end());

    const int idx = static_cast<int>(std::distance(vec.begin(), it));
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
        for (auto& o : offset) {
            o *= -1;
        }
    }
    return offset;
}

#ifdef AMR_WIND_USE_NETCDF
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
#endif

} // namespace

void InletData::resize(const int size)
{
    m_data_n.resize(size);
    m_data_np1.resize(size);
    m_data_interp.resize(size);
}

void InletData::define_plane(const amrex::Orientation ori)
{
    m_data_n[ori] = std::make_unique<PlaneVector>();
    m_data_np1[ori] = std::make_unique<PlaneVector>();
    m_data_interp[ori] = std::make_unique<PlaneVector>();
}

void InletData::define_level_data(
    const amrex::Orientation ori, const amrex::Box& bx, const size_t nc)
{
    if (!this->is_populated(ori)) {
        return;
    }
    m_data_n[ori]->push_back(amrex::FArrayBox(bx, static_cast<int>(nc)));
    m_data_np1[ori]->push_back(amrex::FArrayBox(bx, static_cast<int>(nc)));
    m_data_interp[ori]->push_back(amrex::FArrayBox(bx, static_cast<int>(nc)));
}

#ifdef AMR_WIND_USE_NETCDF
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
    AMREX_ALWAYS_ASSERT(((m_tn <= time) && (time <= m_tnp1)));

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

#endif

void InletData::read_data_native(
    const amrex::OrientationIter oit,
    amrex::BndryRegister& bndry_n,
    amrex::BndryRegister& bndry_np1,
    const int lev,
    const Field* fld,
    const amrex::Real time,
    const amrex::Vector<amrex::Real>& times)
{
    const size_t nc = fld->num_comp();
    const int nstart =
        static_cast<int>(m_components[static_cast<int>(fld->id())]);

    const int idx = closest_index(times, time);
    const int idxp1 = idx + 1;

    m_tn = times[idx];
    m_tnp1 = times[idxp1];

    auto ori = oit();

    AMREX_ALWAYS_ASSERT(((m_tn <= time) && (time <= m_tnp1)));
    AMREX_ALWAYS_ASSERT(fld->num_comp() == bndry_n[ori].nComp());
    AMREX_ASSERT(bndry_n[ori].boxArray() == bndry_np1[ori].boxArray());

    const int normal = ori.coordDir();
    const auto& bbx = (*m_data_n[ori])[lev].box();
    const amrex::IntVect v_offset = offset(ori.faceDir(), normal);

    amrex::MultiFab bndry(
        bndry_n[ori].boxArray(), bndry_n[ori].DistributionMap(),
        bndry_n[ori].nComp(), 0, amrex::MFInfo());

    for (amrex::MFIter mfi(bndry); mfi.isValid(); ++mfi) {

        const auto& vbx = mfi.validbox();
        const auto& bndry_n_arr = bndry_n[ori].array(mfi);
        const auto& bndry_arr = bndry.array(mfi);

        const auto& bx = bbx & vbx;
        if (bx.isEmpty()) {
            continue;
        }

        amrex::ParallelFor(
            bx, nc, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                bndry_arr(i, j, k, n) =
                    0.5 *
                    (bndry_n_arr(i, j, k, n) +
                     bndry_n_arr(
                         i + v_offset[0], j + v_offset[1], k + v_offset[2], n));
            });
    }

    bndry.copyTo((*m_data_n[ori])[lev], 0, nstart, static_cast<int>(nc));

    for (amrex::MFIter mfi(bndry); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto& bndry_np1_arr = bndry_np1[ori].array(mfi);
        const auto& bndry_arr = bndry.array(mfi);

        const auto& bx = bbx & vbx;
        if (bx.isEmpty()) {
            continue;
        }

        amrex::ParallelFor(
            bx, nc, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                bndry_arr(i, j, k, n) =
                    0.5 *
                    (bndry_np1_arr(i, j, k, n) +
                     bndry_np1_arr(
                         i + v_offset[0], j + v_offset[1], k + v_offset[2], n));
            });
    }

    bndry.copyTo((*m_data_np1[ori])[lev], 0, nstart, static_cast<int>(nc));
}

void InletData::interpolate(const amrex::Real time)
{
    m_tinterp = time;
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        if (!this->is_populated(ori)) {
            continue;
        }

        const int lnlevels = static_cast<int>(m_data_n[ori]->size());
        for (int lev = 0; lev < lnlevels; ++lev) {

            const auto& datn = (*m_data_n[ori])[lev];
            const auto& datnp1 = (*m_data_np1[ori])[lev];
            auto& dati = (*m_data_interp[ori])[lev];
            dati.linInterp<amrex::RunOn::Device>(
                datn, 0, datnp1, 0, m_tn, m_tnp1, m_tinterp, datn.box(), 0,
                dati.nComp());
        }
    }
}

bool InletData::is_populated(amrex::Orientation ori) const
{
    return m_data_n[ori] != nullptr;
}

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

    pp.query("bndry_write_frequency", m_write_frequency);
    pp.queryarr("bndry_planes", m_planes);
    pp.query("bndry_output_start_time", m_out_start_time);
    pp.queryarr("bndry_var_names", m_var_names);
    pp.get("bndry_file", m_filename);
    pp.query("bndry_output_format", m_out_fmt);

#ifndef AMR_WIND_USE_NETCDF
    if (m_out_fmt == "netcdf") {
        amrex::Print()
            << "Warning: boundary output format using netcdf must link netcdf "
               "library, changing output to native format"
            << std::endl;
        m_out_fmt = "native";
    }
#endif

    if (!(m_out_fmt == "native" || m_out_fmt == "netcdf")) {
        amrex::Print() << "Warning: boundary output format not recognized, "
                          "changing to native format"
                       << std::endl;
        m_out_fmt = "native";
    }

    // only used for native format
    m_time_file = m_filename + "/time.dat";
}

void ABLBoundaryPlane::post_init_actions()
{
    if (!m_is_initialized) {
        return;
    }
    initialize_data();
    write_header();
    write_file();
    read_header();
    read_file();
}

void ABLBoundaryPlane::pre_advance_work()
{
    if (!m_is_initialized) {
        return;
    }
    read_file();
}

void ABLBoundaryPlane::post_advance_work()
{
    if (!m_is_initialized) {
        return;
    }
    write_file();
}

void ABLBoundaryPlane::initialize_data()
{
    BL_PROFILE("amr-wind::ABLBoundaryPlane::initialize_data");
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
}

void ABLBoundaryPlane::write_header()
{
    BL_PROFILE("amr-wind::ABLBoundaryPlane::write_header");
    if (m_io_mode != io_mode::output) {
        return;
    }

#ifdef AMR_WIND_USE_NETCDF

    if (m_out_fmt == "netcdf") {
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

            const int nlevels = m_repo.num_active_levels();
            for (int lev = 0; lev < nlevels; ++lev) {

                // Only do this if the output plane intersects with data on this
                // level
                const amrex::Box& minBox = m_mesh.boxArray(lev).minimalBox();
                if (!box_intersects_boundary(minBox, lev, ori)) break;

                auto lev_grp = plane_grp.def_group(level_name(lev));
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
                            name, NC_DOUBLE,
                            {"nt", dirs[perp[0]], dirs[perp[1]]});
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

            const int nlevels = plane_grp.num_groups();
            for (int lev = 0; lev < nlevels; ++lev) {
                auto lev_grp = plane_grp.group(level_name(lev));

                const auto& dx = m_mesh.Geom(lev).CellSizeArray();
                const amrex::Box& minBox = m_mesh.boxArray(lev).minimalBox();
                const auto& lo = minBox.loVect();
                const auto& hi = minBox.hiVect();
                const amrex::Vector<amrex::Real> pdx{
                    {dx[perp[0]], dx[perp[1]]}};
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
    }

#endif

    if (amrex::ParallelDescriptor::IOProcessor() && m_out_fmt == "native") {
        // generate time file
        amrex::UtilCreateCleanDirectory(m_filename, false);
        std::ofstream oftime(m_time_file, std::ios::out);
        oftime.close();
    }
}

void ABLBoundaryPlane::write_file()
{
    BL_PROFILE("amr-wind::ABLBoundaryPlane::write_file");
    const amrex::Real time = m_time.new_time();
    const int t_step = m_time.time_index();

    // Only output data if at the desired timestep
    if ((t_step % m_write_frequency != 0) || ((m_io_mode != io_mode::output)) ||
        (time < m_out_start_time)) {
        return;
    }

    for (auto* fld : m_fields) {
        fld->fillpatch(m_time.current_time());
    }

#ifdef AMR_WIND_USE_NETCDF

    if (m_out_fmt == "netcdf") {
        amrex::Print() << "\nWriting NetCDF file " << m_filename << " at time "
                       << time << std::endl;

        auto ncf = ncutils::NCFile::open_par(
            m_filename, NC_WRITE | NC_NETCDF4 | NC_MPIIO,
            amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

        auto v_time = ncf.var("time");
        v_time.par_access(NC_COLLECTIVE);
        const size_t nt = ncf.dim("nt").len();
        v_time.put(&time, {nt}, {1});

        for (amrex::OrientationIter oit; oit; ++oit) {
            auto ori = oit();
            const std::string plane = m_plane_names[ori];

            if (std::find(m_planes.begin(), m_planes.end(), plane) ==
                m_planes.end())
                continue;

            const int nlevels = ncf.group(plane).num_groups();
            for (auto* fld : m_fields) {
                for (int lev = 0; lev < nlevels; ++lev) {
                    auto grp = ncf.group(plane).group(level_name(lev));
                    write_data(grp, ori, lev, fld);
                }
            }
        }

        m_out_counter++;
    }

#endif

    if (m_out_fmt == "native") {
        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ofstream oftime(m_time_file, std::ios::out | std::ios::app);
            oftime << t_step << ' ' << time << '\n';
            oftime.close();
        }

        const std::string chkname =
            m_filename + amrex::Concatenate("/bndry_output", t_step);

        amrex::Print() << "Writing abl boundary checkpoint file " << chkname
                       << " at time " << time << std::endl;

        const std::string level_prefix = "Level_";
        amrex::PreBuildDirectorHierarchy(chkname, level_prefix, 1, true);

        // for now only output level 0
        const int lev = 0;

        for (auto* fld : m_fields) {

            auto& field = *fld;

            const auto& geom = field.repo().mesh().Geom();

            // note: by using the entire domain box we end up using 1 processor
            // to hold all boundaries
            amrex::Box domain = geom[lev].Domain();
            amrex::BoxArray ba(domain);
            amrex::DistributionMapping dm{ba};

            amrex::BndryRegister bndry(
                ba, dm, m_in_rad, m_out_rad, m_extent_rad, field.num_comp());

            bndry.copyFrom(
                field(lev), 0, 0, 0, field.num_comp(), geom[lev].periodicity());

            std::string filename = amrex::MultiFabFileFullPrefix(
                lev, chkname, level_prefix, field.name());

            // print individual faces
            for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
                auto ori = oit();
                const std::string plane = m_plane_names[ori];

                if (std::find(m_planes.begin(), m_planes.end(), plane) ==
                    m_planes.end()) {
                    continue;
                }

                std::string facename =
                    amrex::Concatenate(filename + '_', ori, 1);
                bndry[ori].write(facename);
            }
        }
    }
}

void ABLBoundaryPlane::read_header()
{
    BL_PROFILE("amr-wind::ABLBoundaryPlane::read_header");
    if (m_io_mode != io_mode::input) {
        return;
    }

    // TODO: overallocate this for now
    m_in_data.resize(2 * AMREX_SPACEDIM);

#ifdef AMR_WIND_USE_NETCDF

    if (m_out_fmt == "netcdf") {
        amrex::Print() << "Reading input NetCDF file: " << m_filename
                       << std::endl;
        auto ncf = ncutils::NCFile::open_par(
            m_filename, NC_NOWRITE | NC_NETCDF4 | NC_MPIIO,
            amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

        // Store the input file times and reset to start at 0
        const size_t nt = ncf.dim("nt").len();
        m_in_times.resize(nt);
        ncf.var("time").get(m_in_times.data());

        // Sanity check the input file time
        AMREX_ALWAYS_ASSERT(m_in_times[0] <= m_time.current_time());

        for (auto& plane_grp : ncf.all_groups()) {
            int normal, face_dir;
            plane_grp.var("normal").get(&normal);
            plane_grp.var("side").get(&face_dir);
            const amrex::GpuArray<int, 2> perp = perpendicular_idx(normal);
            const amrex::Orientation ori(
                normal, amrex::Orientation::Side(face_dir));

            m_in_data.define_plane(ori);

            const int nlevels = plane_grp.num_groups();
            for (int lev = 0; lev < nlevels; ++lev) {
                auto lev_grp = plane_grp.group(level_name(lev));

                // sanity checks to ensure grid-to-grid matching
                const amrex::Box& minBox = m_mesh.boxArray(lev).minimalBox();
                const auto& lo = minBox.loVect();
                const auto& hi = minBox.hiVect();
                const auto& dx = m_mesh.Geom(lev).CellSizeArray();
                const amrex::Vector<amrex::Real> pdx{
                    {dx[perp[0]], dx[perp[1]]}};
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
                plo[normal] = ori.isHigh() ? hi[normal] + 1 : -1;
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
    }
#endif

    if (m_out_fmt == "native") {

        int time_file_length = 0;

        if (amrex::ParallelDescriptor::IOProcessor()) {

            std::string line;
            std::ifstream time_file(m_time_file);
            if (!time_file.good()) {
                amrex::Abort("Cannot find time file: " + m_time_file);
            }
            while (std::getline(time_file, line)) {
                ++time_file_length;
            }

            time_file.close();
        }

        amrex::ParallelDescriptor::Bcast(
            &time_file_length, 1,
            amrex::ParallelDescriptor::IOProcessorNumber(),
            amrex::ParallelDescriptor::Communicator());

        m_in_times.resize(time_file_length);
        m_in_timesteps.resize(time_file_length);

        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ifstream time_file(m_time_file);
            for (int i = 0; i < time_file_length; ++i) {
                time_file >> m_in_timesteps[i] >> m_in_times[i];
            }
            time_file.close();
        }

        amrex::ParallelDescriptor::Bcast(
            m_in_timesteps.data(), time_file_length,
            amrex::ParallelDescriptor::IOProcessorNumber(),
            amrex::ParallelDescriptor::Communicator());

        amrex::ParallelDescriptor::Bcast(
            m_in_times.data(), time_file_length,
            amrex::ParallelDescriptor::IOProcessorNumber(),
            amrex::ParallelDescriptor::Communicator());

        int nc = 0;
        for (auto* fld : m_fields) {
            m_in_data.component(static_cast<int>(fld->id())) = nc;
            nc += fld->num_comp();
        }

        // TODO: need to generalize to lev > 0 somehow
        const int lev = 0;
        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            auto ori = oit();

            // TODO: would be safer and less storage to not allocate all of
            // these but we do not use m_planes for input and need to detect
            // mass inflow from field bcs same for define level data below
            m_in_data.define_plane(ori);

            const amrex::Box& minBox = m_mesh.boxArray(lev).minimalBox();

            amrex::IntVect plo(minBox.loVect());
            amrex::IntVect phi(minBox.hiVect());
            const int normal = ori.coordDir();
            plo[normal] = ori.isHigh() ? minBox.hiVect()[normal] + 1 : -1;
            phi[normal] = ori.isHigh() ? minBox.hiVect()[normal] + 1 : -1;
            const amrex::Box pbx(plo, phi);
            m_in_data.define_level_data(ori, pbx, nc);
        }
    }
}

void ABLBoundaryPlane::read_file()
{
    BL_PROFILE("amr-wind::ABLBoundaryPlane::read_file");
    if (m_io_mode != io_mode::input) {
        return;
    }

    // populate planes and interpolate
    const amrex::Real time = m_time.new_time();
    AMREX_ALWAYS_ASSERT((m_in_times[0] <= time) && (time < m_in_times.back()));

    // return early if current data files can still be interpolated in time
    if ((m_in_data.tn() <= time) && (time < m_in_data.tnp1())) {
        m_in_data.interpolate(time);
        return;
    }

#ifdef AMR_WIND_USE_NETCDF
    if (m_out_fmt == "netcdf") {

        auto ncf = ncutils::NCFile::open_par(
            m_filename, NC_NOWRITE | NC_NETCDF4 | NC_MPIIO,
            amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

        for (amrex::OrientationIter oit; oit; ++oit) {
            auto ori = oit();
            if (!m_in_data.is_populated(ori)) continue;

            const std::string plane = m_plane_names[ori];
            const int nlevels = ncf.group(plane).num_groups();
            for (auto* fld : m_fields) {
                for (int lev = 0; lev < nlevels; ++lev) {
                    auto grp = ncf.group(plane).group(level_name(lev));
                    m_in_data.read_data(grp, ori, lev, fld, time, m_in_times);
                }
            }
        }
    }

#endif

    if (m_out_fmt == "native") {

        const int index = closest_index(m_in_times, time);
        const int t_step1 = m_in_timesteps[index];
        const int t_step2 = m_in_timesteps[index + 1];

        AMREX_ALWAYS_ASSERT(
            (m_in_times[index] <= time) && (time <= m_in_times[index + 1]));

        const std::string chkname1 =
            m_filename + amrex::Concatenate("/bndry_output", t_step1);
        const std::string chkname2 =
            m_filename + amrex::Concatenate("/bndry_output", t_step2);

        const std::string level_prefix = "Level_";

        const int lev = 0;
        for (auto* fld : m_fields) {

            auto& field = *fld;
            const auto& geom = field.repo().mesh().Geom();

            amrex::Box domain = geom[lev].Domain();
            amrex::BoxArray ba(domain);
            amrex::DistributionMapping dm{ba};

            amrex::BndryRegister bndry1(
                ba, dm, m_in_rad, m_out_rad, m_extent_rad, field.num_comp());
            amrex::BndryRegister bndry2(
                ba, dm, m_in_rad, m_out_rad, m_extent_rad, field.num_comp());

            bndry1.setVal(1.0e13);
            bndry2.setVal(1.0e13);

            std::string filename1 = amrex::MultiFabFileFullPrefix(
                lev, chkname1, level_prefix, field.name());
            std::string filename2 = amrex::MultiFabFileFullPrefix(
                lev, chkname2, level_prefix, field.name());

            for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
                auto ori = oit();

                if ((!m_in_data.is_populated(ori)) ||
                    (field.bc_type()[ori] != BC::mass_inflow)) {
                    continue;
                }

                std::string facename1 =
                    amrex::Concatenate(filename1 + '_', ori, 1);
                std::string facename2 =
                    amrex::Concatenate(filename2 + '_', ori, 1);

                bndry1[ori].read(facename1);
                bndry2[ori].read(facename2);

                m_in_data.read_data_native(
                    oit, bndry1, bndry2, lev, fld, time, m_in_times);
            }
        }
    }

    m_in_data.interpolate(time);
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void ABLBoundaryPlane::populate_data(
    const int lev,
    const amrex::Real time,
    Field& fld,
    amrex::MultiFab& mfab,
    const int dcomp,
    const int orig_comp) const
{

    BL_PROFILE("amr-wind::ABLBoundaryPlane::populate_data");

    if (m_io_mode != io_mode::input) {
        return;
    }

    AMREX_ALWAYS_ASSERT(
        ((m_in_data.tn() <= time) || (time <= m_in_data.tnp1())));
    AMREX_ALWAYS_ASSERT(std::abs(time - m_in_data.tinterp()) < 1e-12);

    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        if ((!m_in_data.is_populated(ori)) ||
            (fld.bc_type()[ori] != BC::mass_inflow)) {
            continue;
        }

        // Only proceed with data population if fine levels touch the boundary
        if (lev > 0) {
            const amrex::Box& minBox = m_mesh.boxArray(lev).minimalBox();
            if (!box_intersects_boundary(minBox, lev, ori)) {
                continue;
            }
        }

        // Ensure inflow data exists at this level
        if (lev >= m_in_data.nlevels(ori)) {
            amrex::Abort("No inflow data at this level.");
        }

        const size_t nc = mfab.nComp();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(mfab, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {

            const auto& sbx = mfi.growntilebox(1);
            const auto& src = m_in_data.interpolate_data(ori, lev);
            const auto& bx = sbx & src.box();
            if (bx.isEmpty()) {
                continue;
            }

            const auto& dest = mfab.array(mfi);
            const auto& src_arr = src.array();
            const int nstart = m_in_data.component(static_cast<int>(fld.id()));
            amrex::ParallelFor(
                bx, nc,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    dest(i, j, k, n + dcomp) =
                        src_arr(i, j, k, n + nstart + orig_comp);
                });
        }
    }

    const auto& geom = fld.repo().mesh().Geom();
    mfab.EnforcePeriodicity(
        0, mfab.nComp(), amrex::IntVect(1), geom[lev].periodicity());
}

#ifdef AMR_WIND_USE_NETCDF
void ABLBoundaryPlane::write_data(
    const ncutils::NCGroup& grp,
    const amrex::Orientation ori,
    const int lev,
    const Field* fld)
{
    BL_PROFILE("amr-wind::ABLBoundaryPlane::write_data");
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
    const auto& dhi = domain.hiVect();

    AMREX_ALWAYS_ASSERT(dlo[0] == 0 && dlo[1] == 0 && dlo[2] == 0);

    grp.var(name).par_access(NC_COLLECTIVE);

    // TODO optimization
    // - move buffer outside this function, probably best as a member
    // - place in object to access as ori/lev/fld
    // - sizing and start/counts should be done only on init and regrid
    const int n_buffers = m_mesh.boxArray(lev).size();
    amrex::Vector<BufferData> buffers(n_buffers);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi((*fld)(lev), false); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.tilebox();
        const auto& blo = bx.loVect();
        const auto& bhi = bx.hiVect();

        if (blo[normal] == dlo[normal] && ori.isLow()) {
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
        } else if (bhi[normal] == dhi[normal] && ori.isHigh()) {
            amrex::IntVect lo(blo);
            amrex::IntVect hi(bhi);
            // shift by one to reuse impl_buffer_field
            lo[normal] = dhi[normal] + 1;
            hi[normal] = dhi[normal] + 1;
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

//! True if box intersects the boundary
bool ABLBoundaryPlane::box_intersects_boundary(
    const amrex::Box& bx, const int lev, const amrex::Orientation ori) const
{
    const amrex::Box& domBox = m_mesh.Geom(lev).Domain();
    const int normal = ori.coordDir();
    amrex::IntVect plo(domBox.loVect());
    amrex::IntVect phi(domBox.hiVect());
    plo[normal] = ori.isHigh() ? domBox.loVect()[normal] : 0;
    phi[normal] = ori.isHigh() ? domBox.hiVect()[normal] : 0;
    const amrex::Box pbx(plo, phi);
    const auto& intersection = bx & pbx;
    return !intersection.isEmpty();
}

} // namespace amr_wind
