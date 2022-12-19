#include "amr-wind/wind_energy/ABLStats.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/DirectionSelector.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/equation_systems/icns/source_terms/ABLForcing.H"
#include "amr-wind/equation_systems/PDEHelpers.H"

#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"

namespace amr_wind {

namespace {
struct TemperatureGradient
{
    // cppcheck-suppress unusedStructMember
    double grad_z;
    // cppcheck-suppress unusedStructMember
    double max_grad_loc;
};
} // namespace

ABLStats::ABLStats(
    CFDSim& sim, const ABLWallFunction& abl_wall_func, const int dir)
    : m_sim(sim)
    , m_abl_wall_func(abl_wall_func)
    , m_temperature(sim.repo().get_field("temperature"))
    , m_mueff(sim.pde_manager().icns().fields().mueff)
    , m_pa_vel(sim, dir)
    , m_pa_temp(m_temperature, sim.time(), dir)
    , m_pa_vel_fine(sim, dir)
    , m_pa_temp_fine(m_temperature, sim.time(), dir)
    , m_pa_mueff(m_mueff, sim.time(), dir)
    , m_pa_tt(m_pa_temp, m_pa_temp)
    , m_pa_tu(m_pa_vel, m_pa_temp)
    , m_pa_uu(m_pa_vel, m_pa_vel)
    , m_pa_uuu(m_pa_vel, m_pa_vel, m_pa_vel)
{}

ABLStats::~ABLStats() = default;

void ABLStats::post_init_actions()
{
    initialize();
    calc_averages();
}

void ABLStats::initialize()
{
    BL_PROFILE("amr-wind::ABLStats::initialize");

    {
        amrex::ParmParse pp("ABL");
        pp.query("stats_output_frequency", m_out_freq);
        pp.query("stats_output_format", m_out_fmt);
        pp.query("normal_direction", m_normal_dir);
        AMREX_ASSERT((0 <= m_normal_dir) && (m_normal_dir < AMREX_SPACEDIM));
        pp.query("kappa", m_kappa);
        amrex::Vector<amrex::Real> gravity{{0.0, 0.0, -9.81}};
        pp.queryarr("gravity", gravity);
        m_gravity = utils::vec_mag(gravity.data());
        pp.get("reference_temperature", m_ref_theta);
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

    if (m_out_fmt == "netcdf") {
        prepare_netcdf_file();
    } else {
        prepare_ascii_file();
    }
}

void ABLStats::calc_averages()
{
    m_pa_vel();
    m_pa_temp();
    m_pa_vel_fine();
    m_pa_temp_fine();
}

//! Calculate sfs stress averages
void ABLStats::calc_sfs_stress_avgs(
    ScratchField& sfs_stress, ScratchField& t_sfs_stress)
{

    BL_PROFILE("amr-wind::ABLStats::calc_sfs_stress_avgs");

    auto& repo = m_sim.repo();

    const auto& m_vel = repo.get_field("velocity");
    auto gradVel = repo.create_scratch_field(9);
    fvm::gradient(*gradVel, m_vel);

    const auto& alphaeff = repo.get_field(pde_impl::mueff_name("temperature"));
    auto gradT = repo.create_scratch_field(3);
    fvm::gradient(*gradT, m_temperature);

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        for (amrex::MFIter mfi(m_mueff(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& mueff_arr = m_mueff(lev).array(mfi);
            const auto& alphaeff_arr = alphaeff(lev).array(mfi);
            const auto& gradVel_arr = (*gradVel)(lev).array(mfi);
            const auto& gradT_arr = (*gradT)(lev).array(mfi);
            const auto& sfs_arr = sfs_stress(lev).array(mfi);
            const auto& t_sfs_arr = t_sfs_stress(lev).array(mfi);

            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    sfs_arr(i, j, k, 0) =
                        -mueff_arr(i, j, k) *
                        (gradVel_arr(i, j, k, 1) + gradVel_arr(i, j, k, 3));
                    sfs_arr(i, j, k, 1) =
                        -mueff_arr(i, j, k) *
                        (gradVel_arr(i, j, k, 2) + gradVel_arr(i, j, k, 6));
                    sfs_arr(i, j, k, 2) =
                        -mueff_arr(i, j, k) *
                        (gradVel_arr(i, j, k, 5) + gradVel_arr(i, j, k, 7));

                    for (int icomp = 0; icomp < AMREX_SPACEDIM; icomp++) {
                        t_sfs_arr(i, j, k, icomp) =
                            -alphaeff_arr(i, j, k) * gradT_arr(i, j, k, icomp);
                    }
                });
        }
    }
}

void ABLStats::post_advance_work()
{
    BL_PROFILE("amr-wind::ABLStats::post_advance_work");

    // Always compute mean velocity/temperature profiles
    calc_averages();

    const auto& time = m_sim.time();
    const int tidx = time.time_index();
    // Skip processing if it is not an output timestep
    if (!(tidx % m_out_freq == 0)) {
        return;
    }

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

    m_pa_tt();
    m_pa_tu();
    m_pa_uu();
    m_pa_uuu();

    process_output();
}

template <typename h1_dir, typename h2_dir>
void ABLStats::compute_zi(const h1_dir& h1Sel, const h2_dir& h2Sel)
{

    auto gradT = (this->m_sim.repo())
                     .create_scratch_field(3, m_temperature.num_grow()[0]);
    fvm::gradient(*gradT, m_temperature);

    // Only compute zi using coarsest level

    amrex::Gpu::DeviceVector<TemperatureGradient> tgrad(
        m_ncells_h1 * m_ncells_h2);
    auto* tgrad_ptr = tgrad.data();
    {
        const int normal_dir = m_normal_dir;
        const size_t ncells_h1 = m_ncells_h1;
        amrex::Real dnval = m_dn;
        for (amrex::MFIter mfi(m_temperature(0)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& gradT_arr = (*gradT)(0).array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    int h1 = h1Sel(i, j, k);
                    int h2 = h2Sel(i, j, k);
                    if (tgrad_ptr[h2 * ncells_h1 + h1].grad_z <
                        gradT_arr(i, j, k, normal_dir)) {
                        tgrad_ptr[h2 * ncells_h1 + h1].grad_z =
                            gradT_arr(i, j, k, normal_dir);
                        tgrad_ptr[h2 * ncells_h1 + h1].max_grad_loc =
                            (k + 0.5) * dnval;
                    }
                });
        }
    }

    amrex::Vector<TemperatureGradient> temp_grad(m_ncells_h1 * m_ncells_h2);
    amrex::Vector<TemperatureGradient> gtemp_grad(m_ncells_h1 * m_ncells_h2);
    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, tgrad.begin(), tgrad.end(),
        temp_grad.begin());
#ifdef AMREX_USE_MPI
    MPI_Reduce(
        &temp_grad[0], &gtemp_grad[0], m_ncells_h1 * m_ncells_h2,
        MPI_2DOUBLE_PRECISION, MPI_MAXLOC,
        amrex::ParallelDescriptor::IOProcessorNumber(),
        amrex::ParallelDescriptor::Communicator());
#endif

    m_zi = 0.0;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        for (size_t i = 0; i < m_ncells_h1 * m_ncells_h2; i++) {
            m_zi += gtemp_grad[i].max_grad_loc;
        }
        m_zi /=
            (static_cast<double>(m_ncells_h1) *
             static_cast<double>(m_ncells_h2));
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
    m_pa_vel_fine.output_line_average_ascii(
        stat_dir + "/plane_average_velocity_fine.txt", time.time_index(),
        time.current_time());
    m_pa_temp_fine.output_line_average_ascii(
        stat_dir + "/plane_average_temperature_fine.txt", time.time_index(),
        time.current_time());
    m_pa_mueff.output_line_average_ascii(
        stat_dir + "/plane_average_velocity_mueff.txt", time.time_index(),
        time.current_time());
    m_pa_tt.output_line_average_ascii(
        stat_dir + "/second_moment_temperature_temperature.txt",
        time.time_index(), time.current_time());
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
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    amrex::RealArray abl_forcing = {{0.0, 0.0, 0.0}};
    if (m_abl_forcing != nullptr) {
        abl_forcing = m_abl_forcing->abl_forcing();
    }

    double wstar = 0.0;
    auto Q = m_abl_wall_func.mo().surf_temp_flux;
    if (Q > 1e-10) {
        wstar = std::cbrt(m_gravity * Q * m_zi / m_ref_theta);
    }
    auto L = m_abl_wall_func.mo().L;

    std::ofstream outfile;
    outfile.precision(4);
    outfile.open(
        m_ascii_file_name.c_str(), std::ios_base::out | std::ios_base::app);
    // clang-format off
    outfile << time.new_time() << ", "
            << Q << ", "
            << m_abl_wall_func.mo().surf_temp << ", "
            << m_abl_wall_func.utau() << ", "
            << wstar << ", "
            << L << ", "
            << m_zi << ", "
            << abl_forcing[0] << ", "
            << abl_forcing[1] << std::endl;
    // clang-format on
    outfile.close();
}

void ABLStats::prepare_ascii_file()
{
    BL_PROFILE("amr-wind::ABLStats::prepare_ascii_file");
    amrex::Print() << "WARNING: ABLStats: ASCII output will impact performance"
                   << std::endl;

    // Only I/O processor handles this file I/O
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    const std::string stat_dir = "post_processing";
    const std::string sname =
        amrex::Concatenate("abl_statistics", m_sim.time().time_index());

    if (!amrex::UtilCreateDirectory(stat_dir, 0755)) {
        amrex::CreateDirectoryFailed(stat_dir);
    }
    m_ascii_file_name = stat_dir + "/" + sname + ".txt";

    std::ofstream outfile;
    outfile.open(m_ascii_file_name.c_str(), std::ios_base::out);
    outfile << "Time,   Q, Tsurf, ustar,   wstar,   L,   zi, abl_forcing_x, "
               "abl_forcing_y"
            << std::endl;
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
    ncf.def_var("Q", NC_DOUBLE, {nt_name});
    ncf.def_var("Tsurf", NC_DOUBLE, {nt_name});
    ncf.def_var("ustar", NC_DOUBLE, {nt_name});
    ncf.def_var("wstar", NC_DOUBLE, {nt_name});
    ncf.def_var("L", NC_DOUBLE, {nt_name});
    ncf.def_var("zi", NC_DOUBLE, {nt_name});
    ncf.def_var("abl_forcing_x", NC_DOUBLE, {nt_name});
    ncf.def_var("abl_forcing_y", NC_DOUBLE, {nt_name});

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
    grp.def_var("theta'theta'_r", NC_DOUBLE, two_dim);
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
    grp.def_var("u'theta'_sfs", NC_DOUBLE, two_dim);
    grp.def_var("v'theta'_sfs", NC_DOUBLE, two_dim);
    grp.def_var("w'theta'_sfs", NC_DOUBLE, two_dim);
    grp.def_var("u'v'_sfs", NC_DOUBLE, two_dim);
    grp.def_var("u'w'_sfs", NC_DOUBLE, two_dim);
    grp.def_var("v'w'_sfs", NC_DOUBLE, two_dim);

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

    // First calculate sfs stress averages
    auto sfs_stress = m_sim.repo().create_scratch_field("sfs_stress", 3);
    auto t_sfs_stress = m_sim.repo().create_scratch_field("tsfs_stress", 3);
    calc_sfs_stress_avgs(*sfs_stress, *t_sfs_stress);
    ScratchFieldPlaneAveraging pa_sfs(*sfs_stress, m_sim.time(), m_normal_dir);
    pa_sfs();
    ScratchFieldPlaneAveraging pa_tsfs(
        *t_sfs_stress, m_sim.time(), m_normal_dir);
    pa_tsfs();

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
        auto Q = m_abl_wall_func.mo().surf_temp_flux;
        ncf.var("Q").put(&Q, {nt}, {1});
        auto Tsurf = m_abl_wall_func.mo().surf_temp;
        ncf.var("Tsurf").put(&Tsurf, {nt}, {1});
        if (Q > 1e-10) wstar = std::cbrt(m_gravity * Q * m_zi / m_ref_theta);
        ncf.var("wstar").put(&wstar, {nt}, {1});
        double L = m_abl_wall_func.mo().L;
        ncf.var("L").put(&L, {nt}, {1});
        ncf.var("zi").put(&m_zi, {nt}, {1});

        amrex::RealArray abl_forcing = {{0.0, 0.0, 0.0}};
        if (m_abl_forcing != nullptr) {
            abl_forcing = m_abl_forcing->abl_forcing();
        }
        ncf.var("abl_forcing_x").put(&abl_forcing[0], {nt}, {1});
        ncf.var("abl_forcing_y").put(&abl_forcing[1], {nt}, {1});

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
            auto var = grp.var("theta'theta'_r");
            m_pa_tt.line_moment(0, l_vec);
            var.put(l_vec.data(), start, count);
        }

        {
            amrex::Vector<std::string> var_names{
                "u'theta'_r", "v'theta'_r", "w'theta'_r"};
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
            amrex::Vector<std::string> var_names{
                "u'u'u'_r", "v'v'v'_r", "w'w'w'_r"};
            amrex::Vector<int> var_comp{0, 13, 26};
            for (int i = 0; i < var_comp.size(); i++) {
                m_pa_uuu.line_moment(var_comp[i], l_vec);
                auto var = grp.var(var_names[i]);
                var.put(l_vec.data(), start, count);
            }
        }

        {
            amrex::Vector<std::string> var_names{
                "u'v'_sfs", "u'w'_sfs", "v'w'_sfs"};
            for (int i = 0; i < AMREX_SPACEDIM; i++) {
                pa_sfs.line_average(i, l_vec);
                auto var = grp.var(var_names[i]);
                var.put(l_vec.data(), start, count);
            }
        }

        {
            amrex::Vector<std::string> var_names{
                "u'theta'_sfs", "v'theta'_sfs", "w'theta'_sfs"};
            for (int i = 0; i < AMREX_SPACEDIM; i++) {
                pa_tsfs.line_average(i, l_vec);
                auto var = grp.var(var_names[i]);
                var.put(l_vec.data(), start, count);
            }
        }
    }
    ncf.close();
#endif
}

} // namespace amr_wind
