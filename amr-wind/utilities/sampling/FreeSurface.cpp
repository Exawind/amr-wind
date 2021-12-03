#include "amr-wind/utilities/sampling/FreeSurface.H"
#include "amr-wind/utilities/io_utils.H"
#include <AMReX_MultiFabUtil.H>
#include "amr-wind/utilities/ncutils/nc_interface.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace free_surface {

FreeSurface::FreeSurface(CFDSim& sim, const std::string& label)
    : m_sim(sim), m_label(label), m_vof(sim.repo().get_field("vof"))
{}

FreeSurface::~FreeSurface() = default;

void FreeSurface::initialize()
{
    BL_PROFILE("amr-wind::FreeSurface::initialize");

    {
        amrex::ParmParse pp(m_label);
        pp.query("output_frequency", m_out_freq);
        pp.query("output_format", m_out_fmt);
        // Load parameters of freesurface sampling
        pp.getarr("num_points", m_npts_dir);
        pp.getarr("start", m_start);
        pp.getarr("end", m_end);
        AMREX_ALWAYS_ASSERT(static_cast<int>(m_start.size()) == AMREX_SPACEDIM);
        AMREX_ALWAYS_ASSERT(static_cast<int>(m_end.size()) == AMREX_SPACEDIM);
        AMREX_ALWAYS_ASSERT(static_cast<int>(m_npts_dir.size()) == 2);
    }

    // Calculate total number of points
    m_npts = m_npts_dir[0] * m_npts_dir[1];

    // Turn parameters into 2D grid
    m_locs.resize(m_npts);
    m_out.resize(m_npts);

    // Get size of spacing
    amrex::Vector<amrex::Real> dx = {0.0, 0.0};
    for (int nd = 0; nd < 2; ++nd) {
        int d = m_griddim[nd];
        dx[nd] = (m_end[d] - m_start[d]) / amrex::max(m_npts_dir[nd] - 1, 1);
    }

    // Store locations
    int idx = 0;
    for (int j = 0; j < m_npts_dir[1]; ++j) {
        for (int i = 0; i < m_npts_dir[0]; ++i) {
            m_locs[idx][m_orient] = m_start[m_orient];
            // Initialize output values to initial position
            m_out[idx] = m_start[m_orient];
            for (int nd = 0; nd < 2; ++nd) {
                int d = m_griddim[nd];
                m_locs[idx][d] =
                    m_start[d] + dx[nd] * (i * (1 - nd) + j * (nd));
            }
            ++idx;
        }
    }

    if (m_out_fmt == "netcdf") prepare_netcdf_file();
}

void FreeSurface::post_advance_work()
{

    BL_PROFILE("amr-wind::FreeSurface::post_advance_work");
    const auto& time = m_sim.time();
    const int tidx = time.time_index();
    // Skip processing if it is not an output timestep
    if (!(tidx % m_out_freq == 0)) return;

    // Zero data in output array
    for (int n = 0; n < m_npts; n++) {
        m_out[n] = 0.0;
    }

    // Sum of interface locations at each point (assumes one interface only)
    const int finest_level = m_vof.repo().num_active_levels() - 1;

    for (int lev = 0; lev <= finest_level; lev++) {

        // Use level_mask to identify smallest volume
        amrex::iMultiFab level_mask;
        if (lev < finest_level) {
            level_mask = makeFineMask(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                m_sim.mesh().boxArray(lev + 1), amrex::IntVect(2), 1, 0);
        } else {
            level_mask.define(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                1, 0, amrex::MFInfo());
            level_mask.setVal(1);
        }

        auto& vof = m_vof(lev);
        const auto& geom = m_sim.mesh().Geom(lev);
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
            geom.CellSizeArray();
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo =
            geom.ProbLoArray();

        // Loop points in 2D grid
        for (int n = 0; n < m_npts; n++) {
            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> loc;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                loc[d] = m_locs[n][d];
            }
            m_out[n] = amrex::max(
                m_out[n],
                amrex::ReduceMax(
                    vof, level_mask, 0,
                    [=] AMREX_GPU_HOST_DEVICE(
                        amrex::Box const& bx,
                        amrex::Array4<amrex::Real const> const& vof_arr,
                        amrex::Array4<int const> const& mask_arr)
                        -> amrex::Real {
                        amrex::Real height_fab = 0.0;

                        amrex::Loop(
                            bx, [=, &height_fab](int i, int j, int k) noexcept {
                                // Initialize height measurement
                                amrex::Real ht = plo[2];
                                // Cell location
                                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> xm;
                                xm[0] = plo[0] + (i + 0.5) * dx[0];
                                xm[1] = plo[1] + (j + 0.5) * dx[1];
                                xm[2] = plo[2] + (k + 0.5) * dx[2];
                                // Check if cell contains 2D grid point:
                                // complicated conditional is to avoid
                                // double-counting and includes exception for lo
                                // boundary
                                if (((plo[0] == loc[0] &&
                                      xm[0] - loc[0] == 0.5 * dx[0]) ||
                                     (xm[0] - loc[0] < 0.5 * dx[0] &&
                                      loc[0] - xm[0] <= 0.5 * dx[0])) &&
                                    ((plo[1] == loc[1] &&
                                      xm[1] - loc[1] == 0.5 * dx[1]) ||
                                     (xm[1] - loc[1] < 0.5 * dx[1] &&
                                      loc[1] - xm[1] <= 0.5 * dx[1]))) {
                                    // Check if cell is obviously multiphase,
                                    // then check if cell might have interface
                                    // at top or bottom
                                    if ((vof_arr(i, j, k) < 1.0 &&
                                         vof_arr(i, j, k) > 0.0) ||
                                        (vof_arr(i, j, k) == 0.0 &&
                                         (vof_arr(i, j, k + 1) == 1.0 ||
                                          vof_arr(i, j, k - 1) == 1.0))) {
                                        // Determine which cell to interpolate
                                        // with
                                        if (amrex::max(
                                                vof_arr(i, j, k),
                                                vof_arr(i, j, k + 1)) > 0.5 &&
                                            amrex::min(
                                                vof_arr(i, j, k),
                                                vof_arr(i, j, k + 1)) <= 0.5) {
                                            // Interpolate positive direction
                                            ht = xm[2] +
                                                 (dx[2]) /
                                                     (vof_arr(i, j, k + 1) -
                                                      vof_arr(i, j, k)) *
                                                     (0.5 - vof_arr(i, j, k));
                                        } else {
                                            // Interpolate negative direction
                                            ht = xm[2] -
                                                 (dx[2]) /
                                                     (vof_arr(i, j, k - 1) -
                                                      vof_arr(i, j, k)) *
                                                     (0.5 - vof_arr(i, j, k));
                                        }
                                    }
                                }
                                // Offset by removing lo and contribute to whole
                                height_fab = amrex::max(
                                    height_fab,
                                    mask_arr(i, j, k) * (ht - plo[2]));
                            });
                        return height_fab;
                    }));
        }
    }

    // Loop points in 2D grid
    for (int n = 0; n < m_npts; n++) {
        amrex::ParallelDescriptor::ReduceRealMax(m_out[n]);
    }
    // Add problo back to heights, making them absolute, not relative
    const auto& plo0 = m_sim.mesh().Geom(0).ProbLoArray();
    for (int n = 0; n < m_npts; n++) {
        m_out[n] += plo0[m_orient];
    }

    process_output();
}

void FreeSurface::process_output()
{
    if (m_out_fmt == "ascii") {
        write_ascii();
    } else if (m_out_fmt == "netcdf") {
        write_netcdf();
    } else {
        amrex::Abort("FreeSurface: Invalid output format encountered");
    }
}

void FreeSurface::write_ascii()
{
    BL_PROFILE("amr-wind::FreeSurface::write_ascii");
    amrex::Print()
        << "WARNING: FreeSurface: ASCII output will impact performance"
        << std::endl;

    const std::string post_dir = "post_processing";
    const std::string sname =
        amrex::Concatenate(m_label, m_sim.time().time_index());

    if (!amrex::UtilCreateDirectory(post_dir, 0755)) {
        amrex::CreateDirectoryFailed(post_dir);
    }
    const std::string fname = post_dir + "/" + sname + ".txt";

    if (amrex::ParallelDescriptor::IOProcessor()) {
        //
        // Have I/O processor open file and write everything.
        //
        std::ofstream File;

        File.open(fname.c_str(), std::ios::out | std::ios::trunc);

        if (!File.good()) amrex::FileOpenFailed(fname);

        // Metadata
        File << m_npts << '\n';
        File << m_npts_dir[0] << ' ' << m_npts_dir[1] << '\n';

        // Points in grid
        for (int n = 0; n < m_npts; ++n) {
            File << m_locs[n][0] << ' ' << m_locs[n][1] << ' ' << m_out[n]
                 << '\n';
        }

        File.flush();

        File.close();

        if (!File.good())
            amrex::Abort("FreeSurface::write_ascii(): problem writing file");
    }
}

void FreeSurface::prepare_netcdf_file()
{
#ifdef AMR_WIND_USE_NETCDF

    const std::string post_dir = "post_processing";
    const std::string sname =
        amrex::Concatenate(m_label, m_sim.time().time_index());
    if (!amrex::UtilCreateDirectory(post_dir, 0755)) {
        amrex::CreateDirectoryFailed(post_dir);
    }
    m_ncfile_name = post_dir + "/" + sname + ".nc";

    // Only I/O processor handles NetCDF generation
    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    auto ncf = ncutils::NCFile::create(m_ncfile_name, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "num_time_steps";
    const std::string ngp_name = "num_grid_points";
    ncf.enter_def_mode();
    ncf.put_attr("title", "AMR-Wind data sampling output");
    ncf.put_attr("version", ioutils::amr_wind_version());
    ncf.put_attr("created_on", ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim(ngp_name, m_npts);
    ncf.def_dim("ndim", AMREX_SPACEDIM);
    ncf.def_var("time", NC_DOUBLE, {nt_name});

    // Metadata related to the 2D grid used to sample
    const std::vector<int> ijk{m_npts_dir[0], m_npts_dir[1]};
    ncf.put_attr("ijk_dims", ijk);
    ncf.put_attr("start", m_start);
    ncf.put_attr("end", m_end);

    // Set up array of data for locations
    ncf.def_var("coordinates", NC_DOUBLE, {nt_name, ngp_name, "ndim"});

    ncf.exit_def_mode();

#else
    amrex::Abort(
        "NetCDF support was not enabled during build time. Please "
        "recompile or "
        "use native format");
#endif
}

void FreeSurface::write_netcdf()
{
#ifdef AMR_WIND_USE_NETCDF

    if (!amrex::ParallelDescriptor::IOProcessor()) return;
    auto ncf = ncutils::NCFile::open(m_ncfile_name, NC_WRITE);
    const std::string nt_name = "num_time_steps";
    // Index of the next timestep
    const size_t nt = ncf.dim(nt_name).len();
    {
        auto time = m_sim.time().new_time();
        ncf.var("time").put(&time, {nt}, {1});
    }

    std::vector<size_t> start{nt, 0};
    std::vector<size_t> count{1, 0};

    // Copy m_out into m_locs for simplicity
    for (int n = 0; n < m_npts; ++n) {
        m_locs[n][m_orient] = m_out[n];
    }

    count[1] = m_npts;
    auto var = ncf.var("coordinates");
    var.put(&m_locs[0][0], start, count);

    ncf.close();
#endif
}

} // namespace free_surface
} // namespace amr_wind
