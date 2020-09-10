#include "amr-wind/wind_energy/ABLStats.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/DirectionSelector.H"

#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"

namespace amr_wind {

namespace {
struct TemperatureGradient
{
    double grad_z;
    double max_grad_loc;
};
} // namespace

ABLStats::ABLStats(CFDSim& sim, const ABLWallFunction& abl_wall_func)
    : m_sim(sim)
    , m_abl_wall_func(abl_wall_func)
    , m_temperature(sim.repo().get_field("temperature"))
    , m_mueff(sim.pde_manager().icns().fields().mueff)
    , m_pa_vel(sim, 2)
    , m_pa_temp(m_temperature, sim.time(), 2)
    , m_pa_mueff(m_mueff, sim.time(), 2)
    , m_pa_tu(m_pa_vel, m_pa_temp)
    , m_pa_uu(m_pa_vel, m_pa_vel)
    , m_pa_uuu(m_pa_vel, m_pa_vel, m_pa_vel)
{}

ABLStats::~ABLStats() = default;

void ABLStats::initialize()
{
    BL_PROFILE("amr-wind::ABLStats::initialize");

    {
        amrex::ParmParse pp("ABL");
        pp.query("stats_output_frequency", m_out_freq);
        pp.query("stats_output_format", m_out_fmt);
        pp.query("normal_direction", m_normal_dir);
        AMREX_ASSERT((0 <= m_normal_dir) && (m_normal_dir < AMREX_SPACEDIM));
    }

    // Get normal direction and associated stuff
    const auto& geom = (this->m_sim.repo()).mesh().Geom()[0];
    amrex::Box const& domain = geom.Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);
    switch (m_normal_dir) {
    case 0:
        m_ncells_h1 = dhi.y - dlo.y + 1;
        m_ncells_h2 = dhi.z - dlo.z + 1;
        break;
    case 1:
        m_ncells_h1 = dhi.x - dlo.x + 1;
        m_ncells_h2 = dhi.z - dlo.z + 1;
        break;
    case 2:
        m_ncells_h1 = dhi.x - dlo.x + 1;
        m_ncells_h2 = dhi.y - dlo.y + 1;
        break;
    }
    m_dn = geom.CellSize()[m_normal_dir];

    if (m_out_fmt == "netcdf")
        prepare_netcdf_file();
    else
        prepare_ascii_file();
}

void ABLStats::calc_averages()
{
    m_pa_vel();
    m_pa_temp();
    m_pa_tu();
    m_pa_uu();
    m_pa_uuu();
}

void ABLStats::post_advance_work()
{
    BL_PROFILE("amr-wind::ABLStats::post_advance_work");
    const auto& time = m_sim.time();
    const int tidx = time.time_index();
    // Skip processing if it is not an output timestep
    if (!(tidx % m_out_freq == 0)) return;

    switch (m_normal_dir) {
    case 0:
        compute_zi(YDir(), ZDir());
        break;
    case 1:
        compute_zi(XDir(), ZDir());
        break;
    case 2:
        compute_zi(XDir(), YDir());
        break;
    }

    process_output();
}

template <typename h1_dir, typename h2_dir>
void ABLStats::compute_zi(const h1_dir& h1Sel, const h2_dir& h2Sel)
{

    auto gradT = (this->m_sim.repo())
                     .create_scratch_field(3, m_temperature.num_grow()[0]);
    fvm::gradient(*gradT, m_temperature);

    // Only compute zi using coarsest level

    amrex::Gpu::DeviceVector<TemperatureGradient> tgrad(m_ncells_h1 * m_ncells_h2);
    auto* tgrad_ptr = tgrad.data();
    {
        const int normal_dir = m_normal_dir;
        const int ncells_h2 = m_ncells_h2;
        amrex::Real dnval = m_dn;
        for (amrex::MFIter mfi(m_temperature(0)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& gradT_arr = (*gradT)(0).array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    int h1 = h1Sel(i, j, k);
                    int h2 = h2Sel(i, j, k);
                    if (tgrad_ptr[h2 * ncells_h2 + h1].grad_z <
                        gradT_arr(i, j, k, normal_dir)) {
                        tgrad_ptr[h2 * ncells_h2 + h1].grad_z =
                            gradT_arr(i, j, k, normal_dir);
                        tgrad_ptr[h2 * ncells_h2 + h1].max_grad_loc =
                            (k + 0.5) * dnval;
                    }
                });
        }
    }

    amrex::Vector<TemperatureGradient> temp_grad(m_ncells_h1 * m_ncells_h2);
    amrex::Vector<TemperatureGradient> gtemp_grad(m_ncells_h1 * m_ncells_h2);
    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, tgrad.begin(), tgrad.end(), temp_grad.begin());
#ifdef AMREX_USE_MPI
    MPI_Reduce(
        &temp_grad[0], &gtemp_grad[0], m_ncells_h1 * m_ncells_h2,
        MPI_2DOUBLE_PRECISION, MPI_MAXLOC,
        amrex::ParallelDescriptor::IOProcessorNumber(),
        amrex::ParallelDescriptor::Communicator());
#endif

    m_zi = 0.0;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        for (int i = 0; i < m_ncells_h1 * m_ncells_h2; i++) {
            m_zi += gtemp_grad[i].max_grad_loc;
        }
        m_zi /= (m_ncells_h1 * m_ncells_h2);
    }
}

void ABLStats::process_output()
{

    if (m_out_fmt == "ascii") {
        write_ascii();
    } else if (m_out_fmt == "netcdf") {
        write_netcdf();
    } else {
        amrex::Abort("ABLStats: Invalid output format encountered");
    }
}

void ABLStats::write_ascii()
{
    BL_PROFILE("amr-wind::ABLStats::write_ascii");

    const std::string stat_dir = "post_processing";
    const auto& time = m_sim.time();
    m_pa_vel.output_line_average_ascii(
        stat_dir + "/plane_average_velocity.txt", time.time_index(),
        time.current_time());
    m_pa_temp.output_line_average_ascii(
        stat_dir + "/plane_average_temperature.txt", time.time_index(),
        time.current_time());
    m_pa_mueff.output_line_average_ascii(
        stat_dir + "/plane_average_velocity_mueff.txt", time.time_index(),
        time.current_time());
    m_pa_tu.output_line_average_ascii(
        stat_dir + "/second_moment_temperature_velocity.txt", time.time_index(),
        time.current_time());
    m_pa_uu.output_line_average_ascii(
        stat_dir + "/second_moment_velocity_velocity.txt", time.time_index(),
        time.current_time());
    m_pa_uuu.output_line_average_ascii(
        stat_dir + "/third_moment_velocity_velocity_velocity.txt",
        time.time_index(), time.current_time());

    // Only I/O processor handles this file I/O
    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    std::ofstream outfile;
    outfile.precision(4);
    outfile.open(
        m_ascii_file_name.c_str(), std::ios_base::out | std::ios_base::app);
    outfile << time.new_time() << ", " << m_abl_wall_func.utau() << ", " << 0.0
            << ", " << 0.0 << ", " << m_zi << std::endl;
    outfile.close();
}

void ABLStats::prepare_ascii_file()
{
    BL_PROFILE("amr-wind::ABLStats::prepare_ascii_file");
    amrex::Print() << "WARNING: ABLStats: ASCII output will impact performance"
                   << std::endl;

    // Only I/O processor handles this file I/O
    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    const std::string stat_dir = "post_processing";
    const std::string sname =
        amrex::Concatenate("abl_statistics", m_sim.time().time_index());

    if (!amrex::UtilCreateDirectory(stat_dir, 0755)) {
        amrex::CreateDirectoryFailed(stat_dir);
    }
    m_ascii_file_name = stat_dir + "/" + sname + ".txt";

    std::ofstream outfile;
    outfile.open(m_ascii_file_name.c_str(), std::ios_base::out);
    outfile << "Time,   ustar,   wstar,   L,   zi" << std::endl;
    outfile.close();
}

void ABLStats::prepare_netcdf_file()
{
#ifdef AMR_WIND_USE_NETCDF

    const std::string stat_dir = "post_processing";
    const std::string sname =
        amrex::Concatenate("abl_statistics", m_sim.time().time_index());
    if (!amrex::UtilCreateDirectory(stat_dir, 0755)) {
        amrex::CreateDirectoryFailed(stat_dir);
    }
    m_ncfile_name = stat_dir + "/" + sname + ".nc";

    // Only I/O processor handles NetCDF generation
    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    auto ncf = ncutils::NCFile::create(m_ncfile_name, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "num_time_steps";
    ncf.enter_def_mode();
    ncf.put_attr("title", "AMR-Wind ABL statistics output");
    ncf.put_attr("version", ioutils::amr_wind_version());
    ncf.put_attr("created_on", ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim("ndim", AMREX_SPACEDIM);

    ncf.def_var("time", NC_DOUBLE, {nt_name});
    ncf.def_var("ustar", NC_DOUBLE, {nt_name});
    ncf.def_var("wstar", NC_DOUBLE, {nt_name});
    ncf.def_var("L", NC_DOUBLE, {nt_name});
    ncf.def_var("zi", NC_DOUBLE, {nt_name});

    auto grp = ncf.def_group("mean_profiles");
    size_t n_levels = m_pa_vel.ncell_line();
    const std::string nlevels_name = "nlevels";
    grp.def_dim("nlevels", n_levels);
    const std::vector<std::string> two_dim{nt_name, nlevels_name};
    grp.def_var("h", NC_DOUBLE, {nlevels_name});
    grp.def_var("u", NC_DOUBLE, two_dim);
    grp.def_var("v", NC_DOUBLE, two_dim);
    grp.def_var("w", NC_DOUBLE, two_dim);
    grp.def_var("hvelmag", NC_DOUBLE, two_dim);
    grp.def_var("theta", NC_DOUBLE, two_dim);
    grp.def_var("mueff", NC_DOUBLE, two_dim);
    grp.def_var("u'theta'_r", NC_DOUBLE, two_dim);
    grp.def_var("v'theta'_r", NC_DOUBLE, two_dim);
    grp.def_var("w'theta'_r", NC_DOUBLE, two_dim);
    grp.def_var("u'u'_r", NC_DOUBLE, two_dim);
    grp.def_var("u'v'_r", NC_DOUBLE, two_dim);
    grp.def_var("u'w'_r", NC_DOUBLE, two_dim);
    grp.def_var("v'v'_r", NC_DOUBLE, two_dim);
    grp.def_var("v'w'_r", NC_DOUBLE, two_dim);
    grp.def_var("w'w'_r", NC_DOUBLE, two_dim);
    grp.def_var("u'u'u'_r", NC_DOUBLE, two_dim);
    grp.def_var("v'v'v'_r", NC_DOUBLE, two_dim);
    grp.def_var("w'w'w'_r", NC_DOUBLE, two_dim);

    ncf.exit_def_mode();

    {
        const std::vector<size_t> start{0};
        std::vector<size_t> count{n_levels};
        auto h = grp.var("h");
        h.put(m_pa_vel.line_centroids().data(), start, count);
    }

#else
    amrex::Abort(
        "NetCDF support was not enabled during build time. Please recompile or "
        "use native format");
#endif
}

void ABLStats::write_netcdf()
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
        auto ustar = m_abl_wall_func.utau();
        ncf.var("ustar").put(&ustar, {nt}, {1});
        double wstar = 0.0;
        ncf.var("wstar").put(&wstar, {nt}, {1});
        double L = 0.0;
        ncf.var("L").put(&L, {nt}, {1});
        ncf.var("zi").put(&m_zi, {nt}, {1});

        auto grp = ncf.group("mean_profiles");
        size_t n_levels = m_pa_vel.ncell_line();
        amrex::Vector<amrex::Real> l_vec(n_levels);
        std::vector<size_t> start{nt, 0};
        std::vector<size_t> count{1, n_levels};

        {
            amrex::Vector<std::string> var_names{"u", "v", "w"};
            for (int i = 0; i < AMREX_SPACEDIM; i++) {
                m_pa_vel.line_average(i, l_vec);
                auto var = grp.var(var_names[i]);
                var.put(l_vec.data(), start, count);
            }
        }

        {
            auto var = grp.var("hvelmag");
            var.put(m_pa_vel.line_hvelmag_average().data(), start, count);
        }

        {
            auto var = grp.var("theta");
            var.put(m_pa_temp.line_average().data(), start, count);
        }

        {
            auto var = grp.var("mueff");
            var.put(m_pa_mueff.line_average().data(), start, count);
        }

        {
            amrex::Vector<std::string> var_names{"u'theta'_r", "v'theta'_r",
                                                 "w'theta'_r"};
            for (int i = 0; i < AMREX_SPACEDIM; i++) {
                m_pa_tu.line_moment(i, l_vec);
                auto var = grp.var(var_names[i]);
                var.put(l_vec.data(), start, count);
            }
        }

        {
            amrex::Vector<std::string> var_names{"u'u'_r", "u'v'_r", "u'w'_r",
                                                 "v'v'_r", "v'w'_r", "w'w'_r"};
            amrex::Vector<int> var_comp{0, 1, 2, 4, 5, 8};
            for (int i = 0; i < var_comp.size(); i++) {
                m_pa_uu.line_moment(var_comp[i], l_vec);
                auto var = grp.var(var_names[i]);
                var.put(l_vec.data(), start, count);
            }
        }

        {
            amrex::Vector<std::string> var_names{"u'u'u'_r", "v'v'v'_r",
                                                 "w'w'w'_r"};
            amrex::Vector<int> var_comp{0, 13, 26};
            for (int i = 0; i < var_comp.size(); i++) {
                m_pa_uuu.line_moment(var_comp[i], l_vec);
                auto var = grp.var(var_names[i]);
                var.put(l_vec.data(), start, count);
            }
        }
    }
    ncf.close();
#endif
}

} // namespace amr_wind
