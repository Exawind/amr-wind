#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/ocean_waves/utils/wave_utils_K.H"
#include "amr-wind/ocean_waves/OceanWaves.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"

namespace amr_wind_tests {

class OceanWavesOpTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{32, 4, 4}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 4);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
            amrex::Vector<amrex::Real> probhi{{10.0, 1.0, 1.0}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
        {
            // Periodicity
            amrex::ParmParse pp("geometry");
            amrex::Vector<int> periodic{{0, 0, 0}};
            pp.addarr("is_periodic", periodic);
            // Boundary conditions
            amrex::ParmParse ppxlo("xlo");
            ppxlo.add("type", (std::string) "slip_wall");
            amrex::ParmParse ppylo("ylo");
            ppylo.add("type", (std::string) "slip_wall");
            amrex::ParmParse ppzlo("zlo");
            ppzlo.add("type", (std::string) "slip_wall");
            amrex::ParmParse ppxhi("xhi");
            ppxhi.add("type", (std::string) "slip_wall");
            amrex::ParmParse ppyhi("yhi");
            ppyhi.add("type", (std::string) "slip_wall");
            amrex::ParmParse ppzhi("zhi");
            ppzhi.add("type", (std::string) "slip_wall");
        }
        {
            // Physics
            amrex::ParmParse pp("incflo");
            amrex::Vector<std::string> phystr{"MultiPhase", "OceanWaves"};
            pp.addarr("physics", phystr);
            pp.add("use_godunov", (int)1);
        }
    }
};

namespace {

// Functions for HOS distribution, choose periodic in x and y
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
eta_def(amrex::Real x, amrex::Real y)
{
    return (0.75 + 0.25 * std::sin(2.0 * M_PI * x) * std::cos(2.0 * M_PI * y));
}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
u_def(amrex::Real x, amrex::Real y, amrex::Real z)
{
    return (0.1 * std::sin(2.0 * M_PI * x) + 0.0 * y + 0.5 * z);
}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
v_def(amrex::Real x, amrex::Real y, amrex::Real z)
{
    return (
        0.05 * std::sin(2.0 * M_PI * x) + 0.6 * std::cos(2.0 * M_PI * y) +
        0.0 * z);
}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
w_def(amrex::Real x, amrex::Real y, amrex::Real z)
{
    return (0.0 * x + 0.3 * std::cos(2.0 * M_PI * y) + 0.1 * z);
}

void initialize_relaxation_zone_field(
    const amrex::Box& bx,
    const amrex::Array4<amrex::Real>& theor_farr,
    amrex::Real dx,
    amrex::Real xlo,
    amrex::Real gen_length)
{
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real x = xlo + (i + 0.5) * dx;
        amrex::Real xtilde = std::max(std::min(1. - x / gen_length, 1.0), 0.0);
        theor_farr(i, j, k) =
            std::expm1(std::pow(xtilde, 3.5)) / std::expm1(1.0);
    });
}

void init_relaxation_field(amr_wind::Field& theor_field, amrex::Real gen_length)
{
    const auto& geom = theor_field.repo().mesh().Geom();
    run_algorithm(theor_field, [&](const int lev, const amrex::MFIter& mfi) {
        auto theor_field_arr = theor_field(lev).array(mfi);
        const auto& bx = mfi.validbox();
        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();
        initialize_relaxation_zone_field(
            bx, theor_field_arr, dx[0], problo[0], gen_length);
    });
}

void apply_relaxation_zone_field(
    amr_wind::Field& comp, amr_wind::Field& targ, amrex::Real gen_length)
{

    const auto& geom = comp.repo().mesh().Geom();

    for (int lev = 0; lev < comp.repo().num_active_levels(); ++lev) {
        for (amrex::MFIter mfi(comp(lev)); mfi.isValid(); ++mfi) {
            const auto& gbx = mfi.growntilebox(2);
            const auto& dx = geom[lev].CellSizeArray();
            const auto& problo = geom[lev].ProbLoArray();
            const auto& probhi = geom[lev].ProbHiArray();

            auto comp_arr = comp(lev).array(mfi);
            auto targ_arr = targ(lev).array(mfi);

            amrex::ParallelFor(
                gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = amrex::min(
                        amrex::max(problo[0] + (i + 0.5) * dx[0], problo[0]),
                        probhi[0]);
                    if (x <= problo[0] + gen_length) {
                        const amrex::Real Gamma =
                            amr_wind::ocean_waves::utils::Gamma_generate(
                                x - problo[0], gen_length);
                        comp_arr(i, j, k) = targ_arr(i, j, k) * (1. - Gamma) +
                                            comp_arr(i, j, k) * Gamma;
                    }
                });
        }
    }
}

amrex::Real field_error(amr_wind::Field& comp, amr_wind::Field& targ, int ncomp)
{
    amrex::Real error_total = 0.0;
    int nc = ncomp;

    for (int lev = 0; lev < comp.repo().num_active_levels(); ++lev) {
        error_total += amrex::ReduceSum(
            comp(lev), targ(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& comp_arr,
                amrex::Array4<amrex::Real const> const& targ_arr)
                -> amrex::Real {
                amrex::Real error = 0.0;

                amrex::Loop(
                    bx, nc, [=, &error](int i, int j, int k, int n) noexcept {
                        error += std::abs(
                            comp_arr(i, j, k, n) - targ_arr(i, j, k, n));
                    });

                return error;
            });
    }
    return error_total;
}

amrex::Real field_error(amr_wind::Field& comp, amr_wind::Field& targ)
{
    return field_error(comp, targ, 1);
}

amrex::Real gas_velocity_error(
    amr_wind::Field& vel, amr_wind::Field& vof, amrex::Real gas_vel)
{
    amrex::Real error_total = 0.0;
    const amrex::Real gvel = gas_vel;
    const int nc = 3;

    for (int lev = 0; lev < vel.repo().num_active_levels(); ++lev) {
        error_total += amrex::ReduceSum(
            vel(lev), vof(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& vel_arr,
                amrex::Array4<amrex::Real const> const& vof_arr)
                -> amrex::Real {
                amrex::Real error = 0.0;

                amrex::Loop(
                    bx, nc, [=, &error](int i, int j, int k, int n) noexcept {
                        error +=
                            (vof_arr(i, j, k) < 1e-12
                                 ? std::abs(vel_arr(i, j, k, n) - gvel)
                                 : 0.0);
                    });

                return error;
            });
    }
    return error_total;
}

void write_HOS_txt(std::string HOS_fname, amrex::Real factor)
{
    std::ofstream os(HOS_fname);
    // Write metadata to file
    os << "HOS Time = 0.0\n";
    os << "HOS dt = 0.5\n";
    os << "nx = 32, Lx = 10.0\n";
    os << "ny = 4, Ly = 1.0\n";
    os << "nz = 4, zmin = 0.0, zmax = 1.0\n";
    os << "key\n";
    // Loop and write eta, u, v, w
    for (int i = 0; i < 32; ++i) {
        for (int j = 0; j < 4; ++j) {
            // Interface height at current i, j
            amrex::Real x = (0.5 + i) / 3.2;
            amrex::Real y = (0.5 + j) / 4.0;
            os << factor * eta_def(x, y) << std::endl;
            // Velocity at i,j,k
            for (int k = 0; k < 4; ++k) {
                amrex::Real z = (0.5 + k) / 4.0;
                os << factor * u_def(x, y, z) << " " << factor * v_def(x, y, z)
                   << " " << factor * w_def(x, y, z) << std::endl;
            }
        }
    }
}

void init_reference_fields(
    amr_wind::Field& ref_lvs, amr_wind::Field& ref_vel, amrex::Real mfactor)
{
    const auto& geom = ref_lvs.repo().mesh().Geom();
    amrex::Real fac = mfactor;
    run_algorithm(ref_lvs, [&](const int lev, const amrex::MFIter& mfi) {
        auto lvs_arr = ref_lvs(lev).array(mfi);
        auto vel_arr = ref_vel(lev).array(mfi);
        const auto& bx = mfi.validbox();
        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
            const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
            amrex::Real z = problo[2] + (k + 0.5) * dx[2];
            const amrex::Real eta = fac * eta_def(x, y);
            lvs_arr(i, j, k) = eta - z;
            if (lvs_arr(i, j, k) + 0.5 * dx[2] < 0.0) {
                vel_arr(i, j, k, 0) = 0.0;
                vel_arr(i, j, k, 1) = 0.0;
                vel_arr(i, j, k, 2) = 0.0;
            } else {
                // adjust z for partially liquid cell
                if (std::abs(lvs_arr(i, j, k)) - 0.5 * dx[2] < 0) {
                    z -= 0.5 * lvs_arr(i, j, k);
                }
                vel_arr(i, j, k, 0) = fac * u_def(x, y, z);
                vel_arr(i, j, k, 1) = fac * v_def(x, y, z);
                vel_arr(i, j, k, 2) = fac * w_def(x, y, z);
            }
        });
    });
}

void interp_reference_fields(
    amr_wind::Field& ref_lvs,
    amr_wind::Field& ref_vel,
    amr_wind::Field& new_lvs,
    amr_wind::Field& new_vel,
    amrex::Real tfactor)
{
    const amrex::Real tf = tfactor;
    run_algorithm(ref_lvs, [&](const int lev, const amrex::MFIter& mfi) {
        auto lvs_arr = ref_lvs(lev).array(mfi);
        auto vel_arr = ref_vel(lev).array(mfi);
        auto lvs_new = new_lvs(lev).array(mfi);
        auto vel_new = new_vel(lev).array(mfi);
        const auto& bx = mfi.validbox();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            lvs_arr(i, j, k) += (lvs_new(i, j, k) - lvs_arr(i, j, k)) * tf;
            vel_arr(i, j, k, 0) +=
                (vel_new(i, j, k, 0) - vel_arr(i, j, k, 0)) * tf;
            vel_arr(i, j, k, 1) +=
                (vel_new(i, j, k, 1) - vel_arr(i, j, k, 1)) * tf;
            vel_arr(i, j, k, 2) +=
                (vel_new(i, j, k, 2) - vel_arr(i, j, k, 2)) * tf;
        });
    });
}

} // namespace

TEST_F(OceanWavesOpTest, relaxation_zone)
{
    constexpr double tol = 1.0e-12;

    populate_parameters();
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<int> periodic{{0, 0, 0}};
        pp.addarr("is_periodic", periodic);
    }

    initialize_mesh();

    auto& repo = sim().repo();
    const int ncomp = 1;
    const int nghost = 3;
    amrex::Real gen_length = 4.0;
    auto& comp_field = repo.declare_field("comp_field", ncomp, nghost);
    auto& target_field = repo.declare_field("target_field", ncomp, nghost);
    auto& theoretical_field =
        repo.declare_field("theoretical_field", ncomp, nghost);
    comp_field.setVal(0.0);
    target_field.setVal(1.0);
    init_relaxation_field(theoretical_field, gen_length);
    apply_relaxation_zone_field(comp_field, target_field, gen_length);
    amrex::Real error_total = field_error(comp_field, theoretical_field);
    EXPECT_NEAR(error_total, 0.0, tol);
}

TEST_F(OceanWavesOpTest, gas_phase)
{
    // Write HOS file
    write_HOS_txt("HOSGridData_lev0_0.txt", 1.0);

    constexpr double tol = 1.0e-3;

    populate_parameters();
    {
        // Ocean Waves details
        amrex::ParmParse pp("OceanWaves");
        pp.add("label", (std::string) "lin_ow");
        amrex::ParmParse ppow("OceanWaves.lin_ow");
        ppow.add("type", (std::string) "LinearWaves");
        ppow.add("wave_height", 0.02);
        ppow.add("wave_length", 1.0);
        ppow.add("water_depth", 0.5);
        // Wave generation and numerical beach
        ppow.add("relax_zone_gen_length", 2.0);
        ppow.add("numerical_beach_length", 4.0);
    }
    {
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", 0.1);
    }

    initialize_mesh();

    // ICNS must be initialized for MultiPhase physics, which is needed for
    // OceanWaves
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    // Initialize physics
    sim().init_physics();
    auto& oceanwaves =
        sim().physics_manager().get<amr_wind::ocean_waves::OceanWaves>();
    // Initialize fields
    oceanwaves.pre_init_actions();
    auto& repo = sim().repo();
    for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
        oceanwaves.initialize_fields(lev, mesh().Geom(lev));
    }

    // Modify velocity field
    auto& velocity = repo.get_field("velocity");
    const amrex::Real gas_vel = 1.0;
    velocity.setVal(gas_vel);

    // Do post-init step, which modifies velocity and vof fields
    oceanwaves.post_init_actions();

    // Get vof field
    auto& vof = repo.get_field("vof");

    // Check velocity field to confirm not modified
    amrex::Real error_total = gas_velocity_error(velocity, vof, gas_vel);
    EXPECT_NEAR(error_total, 0.0, tol);
}

TEST_F(OceanWavesOpTest, HOS_init)
{
    // Write HOS file
    write_HOS_txt("HOSGridData_lev0_0.txt", 1.0);

    constexpr double tol = 1.0e-3;

    populate_parameters();
    {
        // Ocean Waves details
        amrex::ParmParse pp("OceanWaves");
        pp.add("label", (std::string) "HOS_ow");
        amrex::ParmParse ppow("OceanWaves.HOS_ow");
        ppow.add("type", (std::string) "HOSWaves");
        ppow.add("HOS_files_prefix", (std::string) "HOSGridData");
        ppow.add("initialize_wave_field", (bool)true);
    }
    {
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", 0.1);
    }

    initialize_mesh();

    // ICNS must be initialized for MultiPhase physics, which is needed for
    // OceanWaves
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    // Initialize physics
    sim().init_physics();
    auto& oceanwaves =
        sim().physics_manager().get<amr_wind::ocean_waves::OceanWaves>();
    // Do initial steps with ocean waves
    oceanwaves.pre_init_actions();
    auto& repo = sim().repo();
    for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
        oceanwaves.initialize_fields(lev, mesh().Geom(lev));
    }

    // Create reference fields
    const int nghost = 3;
    auto& ref_levelset = repo.declare_field("ref_levelset", 1, nghost);
    auto& ref_velocity = repo.declare_field("ref_velocity", 3, nghost);
    auto& ref2_levelset = repo.declare_field("ref2_levelset", 1, nghost);
    auto& ref2_velocity = repo.declare_field("ref2_velocity", 3, nghost);
    // Initialize reference fields
    init_reference_fields(ref_levelset, ref_velocity, 1.0);

    // Check physical initialized fields
    auto& levelset = repo.get_field("levelset");
    auto& velocity = repo.get_field("velocity");
    amrex::Real error_total = field_error(levelset, ref_levelset);
    EXPECT_NEAR(error_total, 0.0, tol);
    error_total = field_error(ref_velocity, velocity, 3);
    EXPECT_NEAR(error_total, 0.0, tol);

    // Do post-init step, which stores initial hos and ow fields
    oceanwaves.post_init_actions();

    // Check HOS fields, should match file at n = 0
    auto& hos_levelset = repo.get_field("hos_levelset");
    auto& hos_velocity = repo.get_field("hos_velocity");
    error_total = field_error(ref_levelset, hos_levelset);
    EXPECT_NEAR(error_total, 0.0, tol);
    error_total = field_error(ref_velocity, hos_velocity, 3);
    EXPECT_NEAR(error_total, 0.0, tol);

    // Check ow fields, should match file at n = 0
    auto& ow_levelset = repo.get_field("ow_levelset");
    auto& ow_velocity = repo.get_field("ow_velocity");
    error_total = field_error(ref_levelset, ow_levelset);
    EXPECT_NEAR(error_total, 0.0, tol);
    error_total = field_error(ref_velocity, ow_velocity, 3);
    EXPECT_NEAR(error_total, 0.0, tol);

    // Create file at n = 1 timestep with slightly different values
    write_HOS_txt("HOSGridData_lev0_1.txt", 0.9);
    // Modify reference fields
    init_reference_fields(ref2_levelset, ref2_velocity, 0.9);

    // Advance time
    sim().time().new_timestep();
    // Update relaxation zones to populate HOS fields
    oceanwaves.post_advance_work();

    // Check HOS fields, should match file at n = 1
    error_total = field_error(ref2_levelset, hos_levelset);
    EXPECT_NEAR(error_total, 0.0, tol);
    error_total = field_error(ref2_velocity, hos_velocity, 3);
    EXPECT_NEAR(error_total, 0.0, tol);

    // Check ow fields, should be interpolated in time
    // dt_sim / dt_HOS = 0.1 / 0.5 = 0.2
    interp_reference_fields(
        ref_levelset, ref_velocity, ref2_levelset, ref2_velocity, 0.2);
    error_total = field_error(ref_levelset, ow_levelset);
    EXPECT_NEAR(error_total, 0.0, tol);
    error_total = field_error(ref_velocity, ow_velocity, 3);
    EXPECT_NEAR(error_total, 0.0, tol);
    // Clean up txt files
    {
        const char* fname = "HOSGridData_lev0_0.txt";
        std::ifstream f(fname);
        if (f.good()) {
            remove(fname);
        }
        // Check that file is removed
        std::ifstream ff(fname);
        EXPECT_FALSE(ff.good());
    }
    {
        const char* fname = "HOSGridData_lev0_1.txt";
        std::ifstream f(fname);
        if (f.good()) {
            remove(fname);
        }
        // Check that file is removed
        std::ifstream ff(fname);
        EXPECT_FALSE(ff.good());
    }
}

TEST_F(OceanWavesOpTest, HOS_restart)
{
    // Write HOS file
    write_HOS_txt("HOSGridData_lev0_0.txt", 1.0);

    constexpr double tol = 1.0e-3;

    populate_parameters();
    {
        // Ocean Waves details
        amrex::ParmParse pp("OceanWaves");
        pp.add("label", (std::string) "HOS_ow");
        amrex::ParmParse ppow("OceanWaves.HOS_ow");
        ppow.add("type", (std::string) "HOSWaves");
        ppow.add("HOS_files_prefix", (std::string) "HOSGridData");
        ppow.add("initialize_wave_field", (bool)true);
    }
    {
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", 0.1);
    }

    initialize_mesh();

    // ICNS must be initialized for MultiPhase physics, which is needed for
    // OceanWaves
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    // Initialize physics
    sim().init_physics();
    auto& oceanwaves =
        sim().physics_manager().get<amr_wind::ocean_waves::OceanWaves>();
    // Do initial steps with ocean waves
    oceanwaves.pre_init_actions();
    auto& repo = sim().repo();
    for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
        oceanwaves.initialize_fields(lev, mesh().Geom(lev));
    }

    // Create reference fields
    const int nghost = 3;
    auto& ref_levelset = repo.declare_field("ref_levelset", 1, nghost);
    auto& ref_velocity = repo.declare_field("ref_velocity", 3, nghost);
    auto& ref2_levelset = repo.declare_field("ref2_levelset", 1, nghost);
    auto& ref2_velocity = repo.declare_field("ref2_velocity", 3, nghost);
    // Initialize reference fields
    init_reference_fields(ref_levelset, ref_velocity, 1.0);

    // Advance time before post_init, representing a restart not at t=0
    sim().time().new_timestep();
    // Do post-init step, which stores initial hos and ow fields
    oceanwaves.post_init_actions();

    // Check HOS fields, should match file at n = 0
    auto& hos_levelset = repo.get_field("hos_levelset");
    auto& hos_velocity = repo.get_field("hos_velocity");
    amrex::Real error_total = field_error(ref_levelset, hos_levelset);
    EXPECT_NEAR(error_total, 0.0, tol);
    error_total = field_error(ref_velocity, hos_velocity, 3);
    EXPECT_NEAR(error_total, 0.0, tol);

    // Check ow fields, should match file at n = 0
    auto& ow_levelset = repo.get_field("ow_levelset");
    auto& ow_velocity = repo.get_field("ow_velocity");
    error_total = field_error(ref_levelset, ow_levelset);
    EXPECT_NEAR(error_total, 0.0, tol);
    error_total = field_error(ref_velocity, ow_velocity, 3);
    EXPECT_NEAR(error_total, 0.0, tol);

    // Create file at n = 1 timestep with slightly different values
    write_HOS_txt("HOSGridData_lev0_1.txt", 0.9);
    // Modify reference fields
    init_reference_fields(ref2_levelset, ref2_velocity, 0.9);

    // Advance time
    sim().time().new_timestep();
    // Update relaxation zones to populate HOS fields
    oceanwaves.post_advance_work();

    // Check HOS fields, should match file at n = 1
    error_total = field_error(ref2_levelset, hos_levelset);
    EXPECT_NEAR(error_total, 0.0, tol);
    error_total = field_error(ref2_velocity, hos_velocity, 3);
    EXPECT_NEAR(error_total, 0.0, tol);

    // Check ow fields, should be interpolated in time
    // 2 * dt_sim / dt_HOS = 0.2 / 0.5 = 0.4
    interp_reference_fields(
        ref_levelset, ref_velocity, ref2_levelset, ref2_velocity, 0.4);
    error_total = field_error(ref_levelset, ow_levelset);
    EXPECT_NEAR(error_total, 0.0, tol);
    error_total = field_error(ref_velocity, ow_velocity, 3);
    EXPECT_NEAR(error_total, 0.0, tol);
    // Clean up txt files
    {
        const char* fname = "HOSGridData_lev0_0.txt";
        std::ifstream f(fname);
        if (f.good()) {
            remove(fname);
        }
        // Check that file is removed
        std::ifstream ff(fname);
        EXPECT_FALSE(ff.good());
    }
    {
        const char* fname = "HOSGridData_lev0_1.txt";
        std::ifstream f(fname);
        if (f.good()) {
            remove(fname);
        }
        // Check that file is removed
        std::ifstream ff(fname);
        EXPECT_FALSE(ff.good());
    }
}

} // namespace amr_wind_tests
