#include "amr-wind/wind_energy/actuator/turbine/external/ExtTurbIface.H"
#include "amr-wind/wind_energy/actuator/turbine/kynema/kynema_types.H"
#include "amr-wind/wind_energy/actuator/turbine/kynema/kynema_wrapper.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_FileSystem.H"

#include <algorithm>
#include <cmath>

namespace exw_kynema {
void build_turbine(
    kynema::interfaces::TurbineInterfaceBuilder& builder,
    const YAML::Node wio,
    const int n_blades,
    const int n_blade_nodes,
    const int n_tower_nodes)
{
    auto& turbine_builder = builder.Turbine();

    for (auto blade : std::views::iota(0U, static_cast<uint>(n_blades))) {
        auto& blade_builder = turbine_builder.Blade(blade)
                                  .SetElementOrder(n_blade_nodes - 1)
                                  .PrescribedRootMotion(false);

        const auto& wio_blade = wio["components"]["blade"];
        const auto& ref_axis = wio_blade["reference_axis"];
        const auto& blade_twist = wio_blade["outer_shape"]["twist"];
        const auto& inertia_matrix =
            wio_blade["structure"]["elastic_properties"]["inertia_matrix"];
        const auto& stiffness_matrix =
            wio_blade["structure"]["elastic_properties"]["stiffness_matrix"];

        const auto axis_grid = ref_axis["x"]["grid"].as<std::vector<double>>();
        const auto x_values = ref_axis["x"]["values"].as<std::vector<double>>();
        const auto y_values = ref_axis["y"]["values"].as<std::vector<double>>();
        const auto z_values = ref_axis["z"]["values"].as<std::vector<double>>();
        for (auto i : std::views::iota(0U, axis_grid.size())) {
            blade_builder.AddRefAxisPoint(
                axis_grid[i], {x_values[i], y_values[i], z_values[i]},
                kynema::interfaces::components::ReferenceAxisOrientation::Z);
        }

        const auto twist_grid = blade_twist["grid"].as<std::vector<double>>();
        const auto twist_values =
            blade_twist["values"].as<std::vector<double>>();
        for (auto i : std::views::iota(0U, twist_grid.size())) {
            blade_builder.AddRefAxisTwist(
                twist_grid[i], twist_values[i] * std::numbers::pi / 180.);
        }

        const auto k_grid = stiffness_matrix["grid"].as<std::vector<double>>();
        const auto m_grid = inertia_matrix["grid"].as<std::vector<double>>();
        const auto mass = inertia_matrix["mass"].as<std::vector<double>>();
        const auto cm_x = inertia_matrix["cm_x"].as<std::vector<double>>();
        const auto cm_y = inertia_matrix["cm_y"].as<std::vector<double>>();
        const auto i_cp = inertia_matrix["i_cp"].as<std::vector<double>>();
        const auto i_edge = inertia_matrix["i_edge"].as<std::vector<double>>();
        const auto i_flap = inertia_matrix["i_flap"].as<std::vector<double>>();
        const auto i_plr = inertia_matrix["i_plr"].as<std::vector<double>>();

        const auto K11 = stiffness_matrix["K11"].as<std::vector<double>>();
        const auto K12 = stiffness_matrix["K12"].as<std::vector<double>>();
        const auto K13 = stiffness_matrix["K13"].as<std::vector<double>>();
        const auto K14 = stiffness_matrix["K14"].as<std::vector<double>>();
        const auto K15 = stiffness_matrix["K15"].as<std::vector<double>>();
        const auto K16 = stiffness_matrix["K16"].as<std::vector<double>>();
        const auto K22 = stiffness_matrix["K22"].as<std::vector<double>>();
        const auto K23 = stiffness_matrix["K23"].as<std::vector<double>>();
        const auto K24 = stiffness_matrix["K24"].as<std::vector<double>>();
        const auto K25 = stiffness_matrix["K25"].as<std::vector<double>>();
        const auto K26 = stiffness_matrix["K26"].as<std::vector<double>>();
        const auto K33 = stiffness_matrix["K33"].as<std::vector<double>>();
        const auto K34 = stiffness_matrix["K34"].as<std::vector<double>>();
        const auto K35 = stiffness_matrix["K35"].as<std::vector<double>>();
        const auto K36 = stiffness_matrix["K36"].as<std::vector<double>>();
        const auto K44 = stiffness_matrix["K44"].as<std::vector<double>>();
        const auto K45 = stiffness_matrix["K45"].as<std::vector<double>>();
        const auto K46 = stiffness_matrix["K46"].as<std::vector<double>>();
        const auto K55 = stiffness_matrix["K55"].as<std::vector<double>>();
        const auto K56 = stiffness_matrix["K56"].as<std::vector<double>>();
        const auto K66 = stiffness_matrix["K66"].as<std::vector<double>>();
        const auto n_sections = k_grid.size();

        for (auto section : std::views::iota(0U, n_sections)) {
            if (abs(m_grid[section] - k_grid[section]) > 1e-8) {
                throw std::runtime_error(
                    "stiffness and mass matrices not on same grid");
            }
            blade_builder.AddSection(
                m_grid[section],
                {{{mass[section], 0., 0., 0., 0.,
                   -mass[section] * cm_y[section]},
                  {0., mass[section], 0., 0., 0.,
                   mass[section] * cm_x[section]},
                  {0., 0., mass[section], mass[section] * cm_y[section],
                   -mass[section] * cm_x[section], 0.},
                  {0., 0., mass[section] * cm_y[section], i_edge[section],
                   -i_cp[section], 0.},
                  {0., 0., -mass[section] * cm_x[section], -i_cp[section],
                   i_flap[section], 0.},
                  {-mass[section] * cm_y[section],
                   mass[section] * cm_x[section], 0., 0., 0., i_plr[section]}}},
                {{
                    {K11[section], K12[section], K13[section], K14[section],
                     K15[section], K16[section]},
                    {K12[section], K22[section], K23[section], K24[section],
                     K25[section], K26[section]},
                    {K13[section], K23[section], K33[section], K34[section],
                     K35[section], K36[section]},
                    {K14[section], K24[section], K34[section], K44[section],
                     K45[section], K46[section]},
                    {K15[section], K25[section], K35[section], K45[section],
                     K55[section], K56[section]},
                    {K16[section], K26[section], K36[section], K46[section],
                     K56[section], K66[section]},
                }},
                kynema::interfaces::components::ReferenceAxisOrientation::Z);
        }
    }

    auto& tower_builder = turbine_builder.Tower()
                              .SetElementOrder(n_tower_nodes - 1)
                              .PrescribedRootMotion(false);

    const auto& wio_tower = wio["components"]["tower"];
    const auto& wio_tower_diameter = wio_tower["outer_shape"]["outer_diameter"];
    const auto tower_diameter_grid =
        wio_tower_diameter["grid"].as<std::vector<double>>();
    const auto tower_diameter_values =
        wio_tower_diameter["values"].as<std::vector<double>>();
    const auto tower_wall_thickness =
        wio_tower["structure"]["layers"][0]["thickness"]["values"]
            .as<std::vector<double>>();
    const auto tower_material_name =
        wio_tower["structure"]["layers"][0]["material"].as<std::string>();
    const auto& tower_ref_axis = wio_tower["reference_axis"];

    const auto axis_grid =
        tower_ref_axis["x"]["grid"].as<std::vector<double>>();
    const auto x_values =
        tower_ref_axis["x"]["values"].as<std::vector<double>>();
    const auto y_values =
        tower_ref_axis["y"]["values"].as<std::vector<double>>();
    const auto z_values =
        tower_ref_axis["z"]["values"].as<std::vector<double>>();
    for (auto i : std::views::iota(0U, axis_grid.size())) {
        tower_builder.AddRefAxisPoint(
            axis_grid[i], {x_values[i], y_values[i], z_values[i]},
            kynema::interfaces::components::ReferenceAxisOrientation::Z);
    }

    const auto tower_base_position = std::array<double, 7>{
        x_values[0], y_values[0], z_values[0], 1., 0., 0., 0.};
    turbine_builder.SetTowerBasePosition(tower_base_position);
    tower_builder.AddRefAxisTwist(0., 0.).AddRefAxisTwist(1., 0.);

    auto tower_material =
        std::find_if(
            std::begin(wio["materials"]), std::end(wio["materials"]),
            [&tower_material_name](const YAML::Node& node) {
                return node["name"].as<std::string>() == tower_material_name;
            })
            ->as<YAML::Node>();

    auto elastic_modulus = tower_material["E"].as<double>();
    auto shear_modulus = tower_material["G"].as<double>();
    auto poisson_ratio = tower_material["nu"].as<double>();
    auto density = tower_material["rho"].as<double>();

    for (auto i : std::views::iota(0U, tower_diameter_grid.size())) {
        const auto section = kynema::beams::GenerateHollowCircleSection(
            tower_diameter_grid[i], elastic_modulus, shear_modulus, density,
            tower_diameter_values[i], tower_wall_thickness[i], poisson_ratio);

        tower_builder.AddSection(
            tower_diameter_grid[i], section.M_star, section.C_star,
            kynema::interfaces::components::ReferenceAxisOrientation::Z);
    }

    const auto& wio_hub = wio["components"]["hub"];
    turbine_builder.SetAzimuthAngle(0.)
        .SetConeAngle(
            wio_hub["cone_angle"].as<double>() * std::numbers::pi / 180.)
        .SetHubDiameter(wio_hub["diameter"].as<double>())
        .SetRotorApexToHub(0.);

    const auto& wio_drivetrain = wio["components"]["drivetrain"];
    turbine_builder
        .SetTowerAxisToRotorApex(
            wio_drivetrain["outer_shape"]["overhang"].as<double>())
        .SetTowerTopToRotorApex(
            wio_drivetrain["outer_shape"]["distance_tt_hub"].as<double>())
        .SetShaftTiltAngle(
            wio_drivetrain["outer_shape"]["uptilt"].as<double>() *
            std::numbers::pi / 180.);

    const auto drivetrain_mass =
        wio_drivetrain["elastic_properties"]["mass"].as<double>();
    const auto drivetrain_inertia =
        wio_drivetrain["elastic_properties"]["inertia_tt"]
            .as<std::vector<double>>();

    turbine_builder.SetYawBearingInertiaMatrix(
        std::array{
            std::array{drivetrain_mass, 0., 0., 0., 0., 0.},
            std::array{0., drivetrain_mass, 0., 0., 0., 0.},
            std::array{0., 0., drivetrain_mass, 0., 0., 0.},
            std::array{
                0., 0., 0., drivetrain_inertia[0], drivetrain_inertia[3],
                drivetrain_inertia[4]},
            std::array{
                0., 0., 0., drivetrain_inertia[3], drivetrain_inertia[1],
                drivetrain_inertia[5]},
            std::array{
                0., 0., 0., drivetrain_inertia[4], drivetrain_inertia[5],
                drivetrain_inertia[2]}});

    // Get hub mass properties from WindIO
    const auto hub_mass = wio_hub["elastic_properties"]["mass"].as<double>();
    const auto hub_inertia =
        wio_hub["elastic_properties"]["inertia"].as<std::vector<double>>();

    // Set the hub inertia matrix in the turbine builder
    turbine_builder.SetHubInertiaMatrix(
        {{{hub_mass, 0., 0., 0., 0., 0.},
          {0., hub_mass, 0., 0., 0., 0.},
          {0., 0., hub_mass, 0., 0., 0.},
          {0., 0., 0., hub_inertia[0], hub_inertia[3], hub_inertia[4]},
          {0., 0., 0., hub_inertia[3], hub_inertia[1], hub_inertia[5]},
          {0., 0., 0., hub_inertia[4], hub_inertia[5], hub_inertia[2]}}});
}

void build_aero(
    kynema::interfaces::TurbineInterfaceBuilder& builder,
    const YAML::Node wio,
    const amrex::Vector<amrex::Real> vel_pt_distr)
{
    auto& aero_builder = builder.Aerodynamics()
                             .EnableAero()
                             .SetNumberOfAirfoils(1UL)
                             .SetAirfoilToBladeMap(std::array{0UL, 0UL, 0UL});

    const auto& airfoil_io = wio["airfoils"];
    auto aero_sections =
        std::vector<kynema::interfaces::components::AerodynamicSection>{};
    auto id = 0UL;
    const auto max_id = airfoil_io.size() - 1UL;
    long id_apos = 0;
    long id_apos_end = id_apos;
    // If there is another aero point before the next spanwise position,
    // duplicate the aero section If there is no aero point in the current aero
    // section, skip it
    for (const auto& af : airfoil_io) {
        const auto s = af["spanwise_position"].as<double>();
        const auto chord = af["chord"].as<double>();
        const auto twist = af["twist"].as<double>() * std::numbers::pi / 180.;
        const auto section_offset_x = af["section_offset_x"].as<double>();
        const auto section_offset_y = af["section_offset_y"].as<double>();
        const auto aerodynamic_center = af["aerodynamic_center"].as<double>();
        auto aoa = af["polars"][0]["re_sets"][0]["cl"]["grid"]
                       .as<std::vector<double>>();
        std::ranges::transform(aoa, std::begin(aoa), [](auto degrees) {
            return degrees * std::numbers::pi / 180.;
        });
        const auto cl = af["polars"][0]["re_sets"][0]["cl"]["values"]
                            .as<std::vector<double>>();
        const auto cd = af["polars"][0]["re_sets"][0]["cd"]["values"]
                            .as<std::vector<double>>();
        const auto cm = af["polars"][0]["re_sets"][0]["cm"]["values"]
                            .as<std::vector<double>>();

        auto point_between = false;
        const bool point_at_end = (id == max_id);
        if (!point_at_end) {
            const auto s_next =
                airfoil_io[id + 1]["spanwise_position"].as<double>();
            point_between =
                (vel_pt_distr[id_apos] > s && vel_pt_distr[id_apos] <= s_next);
            if (point_between) {
                // Check for more points in between
                id_apos_end = id_apos;
                for (int iae = id_apos + 1; iae < vel_pt_distr.size(); ++iae) {
                    if (vel_pt_distr[iae] <= s_next) {
                        id_apos_end = iae;
                    } else {
                        break;
                    }
                }
            }
        } else {
            // If at end, all remaining points belong to this aerodynamic
            // section
            id_apos_end = vel_pt_distr.size() - 1;
        }

        if (point_between || point_at_end) {
            for (id_apos; id_apos <= id_apos_end; ++id_apos) {
                aero_sections.emplace_back(
                    id_apos, s, chord, section_offset_x, section_offset_y,
                    aerodynamic_center, twist, aoa, cl, cd, cm);
            }
        }
        ++id;
    }

    aero_builder.SetAirfoilSections(0UL, aero_sections);
}
} // namespace exw_kynema

namespace ext_turb {

template <>
ExtTurbIface<KynemaTurbine, KynemaSolverData>::~ExtTurbIface()
{
    //! Do deallocation if necessary
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::parse_inputs(
    const amr_wind::CFDSim& sim, const std::string& inp_name)
{
    amrex::ParmParse pp(inp_name);

    const auto& time = sim.time();
    //! Check that the user has not enabled adaptive timestepping
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !time.adaptive_timestep(),
        "Adaptive time-stepping not supported when Kynema is enabled");

    m_dt_cfd = time.delta_t();

    //! Put input-reading steps here
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::allocate_ext_turbines()
{
    BL_PROFILE("amr-wind::KynemaIface::allocate_turbines");
    int nturbines = static_cast<int>(m_turbine_data.size());
    //!! Do things need to be allocated separately in Kynema ??

    m_is_initialized = true;
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::init_solution(
    const int local_id)
{
    BL_PROFILE("amr-wind::KynemaIface::init_solution");
    AMREX_ALWAYS_ASSERT(local_id < static_cast<int>(m_turbine_data.size()));
    AMREX_ALWAYS_ASSERT(m_is_initialized);

    auto& fi = *m_turbine_data[local_id];
    //!! Set up initial solution
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::get_hub_stats(
    const int local_id)
{
    BL_PROFILE("amr-wind::KynemaIface::get_hub_stats");

    auto& fi = *m_turbine_data[local_id];
}

#ifdef AMR_WIND_USE_KYNEMA

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::prepare_netcdf_file(
    KynemaTurbine& fi)
{
#ifdef AMR_WIND_USE_NETCDF
    BL_PROFILE("amr-wind::KynemaIface::prepare_netcdf_file");
    if (!amrex::UtilCreateDirectory(m_solver_data.m_output_dir, 0755)) {
        amrex::CreateDirectoryFailed(m_solver_data.m_output_dir);
    }

    const std::string fname =
        m_solver_data.m_output_dir + "/" + fi.tlabel + ".nc";

    // Don't overwrite existing
    if (amrex::FileSystem::Exists(fname)) {
        return;
    }

    auto ncf = ncutils::NCFile::create(fname, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "num_time_steps";
    const std::string np_name = "num_vel_points";
    ncf.enter_def_mode();
    ncf.put_attr("title", "AMR-Wind Kynema velocity data");
    ncf.put_attr("AMR-Wind_version", amr_wind::ioutils::amr_wind_version());
    ncf.put_attr("created_on", amr_wind::ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim(np_name, fi.from_cfd.u_Len);
    ncf.def_dim("ndim", AMREX_SPACEDIM);
    ncf.def_var("time", NC_FLOAT, {nt_name});
    ncf.def_var("xco", NC_FLOAT, {np_name});
    ncf.def_var("yco", NC_FLOAT, {np_name});
    ncf.def_var("zco", NC_FLOAT, {np_name});
    ncf.def_var("uvel", NC_FLOAT, {nt_name, np_name});
    ncf.def_var("vvel", NC_FLOAT, {nt_name, np_name});
    ncf.def_var("wvel", NC_FLOAT, {nt_name, np_name});
    ncf.exit_def_mode();

    {
        const auto npts = static_cast<size_t>(fi.from_cfd.u_Len);
        auto xco = ncf.var("xco");
        xco.put(fi.position_at_vel(0), {0}, {npts});
        auto yco = ncf.var("yco");
        yco.put(fi.position_at_vel(1), {0}, {npts});
        auto zco = ncf.var("zco");
        zco.put(fi.position_at_vel(2), {0}, {npts});
    }
#else
    amrex::ignore_unused(fi);
    amrex::OutStream() << "WARNING: KynemaIface: NetCDF support was not "
                          "enabled during compile "
                          "time. KynemaIface cannot support restart."
                       << std::endl;
#endif
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::write_velocity_data(
    const KynemaTurbine& fi)
{
#ifdef AMR_WIND_USE_NETCDF
    BL_PROFILE("amr-wind::KynemaIface::write_velocity_data");
    const std::string fname =
        m_solver_data.m_output_dir + "/" + fi.tlabel + ".nc";
    auto ncf = ncutils::NCFile::open(fname, NC_WRITE);
    const std::string nt_name = "num_time_steps";
    const size_t nt = ncf.dim(nt_name).len();
    const auto npts = static_cast<size_t>(fi.from_cfd.u_Len);

    const double time = fi.time_index * fi.dt_ext;
    ncf.var("time").put(&time, {nt}, {1});
    const auto& uu = fi.from_cfd;
    ncf.var("uvel").put(uu.u, {nt, 0}, {1, npts});
    ncf.var("vvel").put(uu.v, {nt, 0}, {1, npts});
    ncf.var("wvel").put(uu.w, {nt, 0}, {1, npts});
#else
    amrex::ignore_unused(fi);
#endif
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::read_velocity_data(
    KynemaTurbine& fi, const ncutils::NCFile& ncf, const size_t tid)
{
#ifdef AMR_WIND_USE_NETCDF
    const auto nt = static_cast<size_t>(tid);
    const auto npts = static_cast<size_t>(fi.from_cfd.u_Len);

    auto& uu = fi.from_cfd;
    ncf.var("uvel").get(uu.u, {nt, 0}, {1, npts});
    ncf.var("vvel").get(uu.v, {nt, 0}, {1, npts});
    ncf.var("wvel").get(uu.w, {nt, 0}, {1, npts});
#else
    amrex::ignore_unused(fi);
    amrex::Abort(
        "KynemaIface::read_velocity_data: AMR-Wind was not compiled with "
        "NetCDF "
        "support");
#endif
}

#else

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::prepare_netcdf_file(
    KynemaTurbine& /*unused*/)
{}
template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::write_velocity_data(
    const KynemaTurbine& /*unused*/)
{}
template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::read_velocity_data(
    KynemaTurbine& /*unused*/,
    const ncutils::NCFile& /*unused*/,
    const size_t /*unused*/)
{}

#endif

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::do_turbine_step(
    KynemaTurbine& fi)
{
    // individual turbine step
    fi.interface->Step();
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::write_turbine_checkpoint(
    int& tid)
{
    // write checkpoint
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::ext_init_turbine(
    KynemaTurbine& fi)
{
    BL_PROFILE("amr-wind::KynemaIface::init_turbine");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        amrex::FileSystem::Exists(fi.input_file),
        "KynemaIface: Cannot find Kynema input file: " + fi.input_file);

    auto builder = kynema::interfaces::TurbineInterfaceBuilder{};
    builder.Solution()
        .EnableDynamicSolve()
        .SetTimeStep(fi.dt_ext)
        .SetDampingFactor(0.0)
        .SetGravity({0., 0., -9.81})
        .SetMaximumNonlinearIterations(6)
        .SetAbsoluteErrorTolerance(1e-6)
        .SetRelativeErrorTolerance(1e-4);

    const YAML::Node wio = YAML::LoadFile(fi.input_file);

    exw_kynema::build_turbine(
        builder, wio, fi.num_blades, fi.num_blade_elem, fi.num_pts_tower);

    // Create distribution of points from base of blade to end, normalized from
    // 0 to 1
    amrex::Vector<amrex::Real> aero_pts_blade{};
    aero_pts_blade.resize(fi.num_pts_blade);
    // Put loop here to create normalized spanwise positions
    const bool uniform_distr{true};
    if (uniform_distr) {
        const amrex::Real dr = 1.0 / static_cast<amrex::Real>(fi.num_pts_blade);
        //!! Is this the way openfast does it? No points exactly at the root or tip
        aero_pts_blade[0] = 0.5 * dr;
        for (int ir = 1; ir < fi.num_pts_blade; ++ir) {
            aero_pts_blade[ir] = (static_cast<amrex::Real>(ir) + 0.5) * dr;
        }
    }

    exw_kynema::build_aero(builder, wio, aero_pts_blade);

    fi.interface = std::make_unique<kynema::interfaces::TurbineInterface>(builder.Build());

    // Determine the number of substeps for Kynema per CFD timestep
    fi.num_substeps = static_cast<int>(std::floor(fi.dt_cfd / fi.dt_ext));

    AMREX_ALWAYS_ASSERT(fi.num_substeps > 0);
    // Check that the time step sizes are consistent and Kynema advances at an
    // integral multiple of CFD timestep
    double dt_err =
        fi.dt_cfd / (static_cast<double>(fi.num_substeps) * fi.dt_ext) - 1.0;
    if (dt_err > 1.0e-12) {
        amrex::Abort(
            "KynemaIFace: Kynema timestep is not an integral "
            "multiple of CFD timestep");
    }
}

// cppcheck-suppress constParameterReference
// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::ext_replay_turbine(
    KynemaTurbine& fi)
{

    // Do we even do this???
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::ext_restart_turbine(
    KynemaTurbine& fi)
{
    BL_PROFILE("amr-wind::KynemaIface::restart_turbine");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        amrex::FileSystem::Exists(fi.checkpoint_file + ".chkp"),
        "KynemaIface: Cannot find Kynema checkpoint file: " +
            fi.checkpoint_file);

    // Determine the number of substeps for Kynema per CFD timestep
    fi.num_substeps = static_cast<int>(std::floor(fi.dt_cfd / fi.dt_ext));

    AMREX_ALWAYS_ASSERT(fi.num_substeps > 0);
    // Check that the time step sizes are consistent and Kynema advances at an
    // integral multiple of CFD timestep
    double dt_err =
        fi.dt_cfd / (static_cast<double>(fi.num_substeps) * fi.dt_ext) - 1.0;
    if (dt_err > 1.0e-4) {
        amrex::Abort(
            "KynemaIFace: Kynema timestep is not an integral "
            "multiple of CFD timestep");
    }
}

template class ExtTurbIface<KynemaTurbine, KynemaSolverData>;

} // namespace ext_turb
