#include "amr-wind/utilities/sampling/FreeSurface.H"
#include "amr-wind/utilities/io_utils.H"
#include <AMReX_MultiFabUtil.H>
#include "amr-wind/utilities/ncutils/nc_interface.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace free_surface {

amrex::Real get_height (
  const int i, const int j, const int k,
  const amrex::GpuArray<amrex::Real, 3>& vof_arr,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& loc,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& problo,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx)
{
    // Initialize height measurement
    amrex::Real height = problo[2];
    // Cell location
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> xm;
    xm[0] = problo[0] + (i + 0.5) * dx[0];
    xm[1] = problo[1] + (j + 0.5) * dx[1];
    xm[2] = problo[2] + (k + 0.5) * dx[2];
    // Check if cell contains 2D grid point: complicated conditional is to avoid
    // double-counting and includes exception for lo boundary
    if (((problo[0] == loc[0] && xm[0] - loc[0] == 0.5 * dx[0]) ||
        (xm[0] - loc[0] < 0.5 * dx[0] && loc[0] - xm[0] <= 0.5 * dx[0])) &&
        ((problo[1] == loc[1] && xm[1] - loc[1] == 0.5 * dx[1]) ||
        (xm[1] - loc[1] < 0.5 * dx[1] && loc[1] - xm[1] <= 0.5 * dx[1]))) {
        // Check if cell is obviously multiphase, then check if cell might have
        // interface at top or bottom
        if ((vof_arr[1] < 1.0 && vof_arr[1] > 0.0) ||
            (vof_arr[1] == 0.0 && (vof_arr[2] == 1.0 || vof_arr[0] == 1.0))) {
            // Determine which cell to interpolate with
            if (amrex::max(vof_arr[1], vof_arr[2]) > 0.5 &&
                amrex::min(vof_arr[1], vof_arr[2]) <= 0.5) {
                // Interpolate positive direction
                height = xm[2] +
                     (dx[2]) / (vof_arr[2] - vof_arr[1]) * (0.5 - vof_arr[1]);
            } else {
                // Interpolate negative direction
                height = xm[2] -
                     (dx[2]) / (vof_arr[0] - vof_arr[1]) * (0.5 - vof_arr[1]);
            }
        }
    }
    height -= problo[2];
    return height;
}

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
                m_locs[idx][d] = m_start[d] + dx[nd] * (i * (1 - nd) + j * (nd));
            }
            ++idx;
        }
    }

    //if (m_out_fmt == "netcdf") prepare_netcdf_file();
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
            m_out[n] += amrex::ReduceSum(
                vof, level_mask, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& vof_arr,
                    amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                    amrex::Real height_fab = 0.0;

                    amrex::Loop(
                        bx, [=, &height_fab](int i, int j, int k) noexcept {
                            amrex::GpuArray<amrex::Real, 3> vof_3 = {
                                vof_arr(i, j, k - 1), vof_arr(i, j, k),
                                vof_arr(i, j, k + 1)};
                            amrex::Real ht =
                                get_height(i, j, k, vof_3, loc, plo, dx);

                            height_fab += mask_arr(i, j, k) * ht;
                        });
                    return height_fab;
                });
        }
    }

    // Loop points in 2D grid
    for (int n = 0; n < m_npts; n++) {
        amrex::ParallelDescriptor::ReduceRealSum(m_out[n]);
    }

    process_output();
}

void FreeSurface::process_output()
{
    // Add problo back to heights
    if (m_out_fmt == "ascii") {
        // write_ascii();
    } else if (m_out_fmt == "netcdf") {
        // write_netcdf();
    } else {
        amrex::Abort("FreeSurface: Invalid output format encountered");
    }
}

    /*void FreeSurface::write_ascii()
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
        m_scontainer->WriteAsciiFile(fname);
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

        auto ncf =
            ncutils::NCFile::create(m_ncfile_name, NC_CLOBBER | NC_NETCDF4);
        const std::string nt_name = "num_time_steps";
        const std::string npart_name = "num_points";
        const std::vector<std::string> two_dim{nt_name, npart_name};
        ncf.enter_def_mode();
        ncf.put_attr("title", "AMR-Wind data sampling output");
        ncf.put_attr("version", ioutils::amr_wind_version());
        ncf.put_attr("created_on", ioutils::timestamp());
        ncf.def_dim(nt_name, NC_UNLIMITED);
        ncf.def_dim("ndim", AMREX_SPACEDIM);
        ncf.def_var("time", NC_DOUBLE, {nt_name});
        // Define groups for each sampler
        for (const auto& obj : m_samplers) {
            auto grp = ncf.def_group(obj->label());

            grp.def_dim(npart_name, obj->num_points());
            obj->define_netcdf_metadata(grp);
            grp.def_var("coordinates", NC_DOUBLE, {npart_name, "ndim"});
            for (const auto& vname : m_var_names)
                grp.def_var(vname, NC_DOUBLE, two_dim);
        }
        ncf.exit_def_mode();

        {
            const std::vector<size_t> start{0, 0};
            std::vector<size_t> count{0, AMREX_SPACEDIM};
            SamplerBase::SampleLocType locs;
            for (const auto& obj : m_samplers) {
                auto grp = ncf.group(obj->label());
                obj->populate_netcdf_metadata(grp);
                obj->sampling_locations(locs);
                auto xyz = grp.var("coordinates");
                count[0] = obj->num_points();
                xyz.put(&locs[0][0], start, count);
            }
        }

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
        std::vector<double> buf(m_total_particles * m_var_names.size(), 0.0);
        m_scontainer->populate_buffer(buf);

        if (!amrex::ParallelDescriptor::IOProcessor()) return;
        auto ncf = ncutils::NCFile::open(m_ncfile_name, NC_WRITE);
        const std::string nt_name = "num_time_steps";
        // Index of the next timestep
        const size_t nt = ncf.dim(nt_name).len();
        {
            auto time = m_sim.time().new_time();
            ncf.var("time").put(&time, {nt}, {1});
        }

        for (const auto& obj : m_samplers) {
            auto grp = ncf.group(obj->label());
            obj->output_netcdf_data(grp, nt);
        }

        std::vector<size_t> start{nt, 0};
        std::vector<size_t> count{1, 0};

        const int nvars = m_var_names.size();
        for (int iv = 0; iv < nvars; ++iv) {
            start[1] = 0;
            count[1] = 0;
            int offset = iv * m_scontainer->num_sampling_particles();
            for (const auto& obj : m_samplers) {
                count[1] = obj->num_points();
                auto grp = ncf.group(obj->label());
                auto var = grp.var(m_var_names[iv]);
                var.put(&buf[offset], start, count);
                offset += count[1];
            }
        }
        ncf.close();
#endif
}*/

} // namespace free_surface
} // namespace amr_wind
