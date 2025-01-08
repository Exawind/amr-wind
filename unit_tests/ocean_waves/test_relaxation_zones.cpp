#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/ocean_waves/utils/wave_utils_K.H"
#include "amr-wind/ocean_waves/OceanWaves.H"
#include "amr-wind/utilities/constants.H"
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
            amrex::Vector<int> ncell{{64, 4, 16}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 4);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0, -0.25, -1}};
            amrex::Vector<amrex::Real> probhi{{8.0, 0.25, 1}};

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
            pp.add("use_godunov", 1);
        }
    }
};

namespace {

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
                            amr_wind::ocean_waves::utils::gamma_generate(
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
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
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
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    return error_total;
}

void make_target_velocity(
    amr_wind::Field& ow_velocity,
    amr_wind::Field& velocity,
    amr_wind::Field& ow_vof)
{
    for (int lev = 0; lev < ow_vof.repo().num_active_levels(); ++lev) {
        for (amrex::MFIter mfi(ow_vof(lev)); mfi.isValid(); ++mfi) {
            const auto& gbx = mfi.growntilebox(2);

            auto ow_vel_arr = ow_velocity(lev).array(mfi);
            auto vel_arr = velocity(lev).const_array(mfi);
            auto ow_vof_arr = ow_vof(lev).const_array(mfi);

            amrex::ParallelFor(
                gbx, velocity.num_comp(),
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    if (ow_vof_arr(i, j, k) <= amr_wind::constants::LOOSE_TOL) {
                        ow_vel_arr(i, j, k, n) = vel_arr(i, j, k, n);
                    }
                });
        }
    }
}

void make_target_density(
    amr_wind::Field& ow_vof, const amrex::Real rho1, const amrex::Real rho2)
{
    for (int lev = 0; lev < ow_vof.repo().num_active_levels(); ++lev) {
        for (amrex::MFIter mfi(ow_vof(lev)); mfi.isValid(); ++mfi) {
            const auto& gbx = mfi.growntilebox(2);
            auto ow_vof_arr = ow_vof(lev).array(mfi);

            amrex::ParallelFor(
                gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    ow_vof_arr(i, j, k) = rho1 * ow_vof_arr(i, j, k) +
                                          rho2 * (1.0 - ow_vof_arr(i, j, k));
                });
        }
    }
}

amrex::Real bdy_error(amr_wind::Field& comp, amr_wind::Field& targ, int ncomp)
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
                        if (i == 0) {
                            error += std::abs(
                                comp_arr(i - 1, j, k, n) -
                                targ_arr(i - 1, j, k, n));
                        }
                    });

                return error;
            });
    }
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    return error_total;
}

amrex::Real bdy_error(amr_wind::Field& comp, amr_wind::Field& targ)
{
    return bdy_error(comp, targ, 1);
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

    constexpr double tol = 1.0e-3;

    populate_parameters();
    {
        // Ocean Waves details
        amrex::ParmParse pp("OceanWaves");
        pp.add("label", (std::string) "lin_ow");
        amrex::ParmParse ppow("OceanWaves.lin_ow");
        ppow.add("type", (std::string) "LinearWaves");
        ppow.add("wave_height", 0.05);
        ppow.add("wave_length", 1.0);
        ppow.add("water_depth", 1.0);
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

TEST_F(OceanWavesOpTest, boundary_fill)
{

    constexpr double tol = 1.0e-3;
    const amrex::Vector<amrex::Real> gas_vel{{1.0, 0.0, 0.0}};

    populate_parameters();
    {
        // Ocean Waves details
        amrex::ParmParse pp("OceanWaves");
        pp.add("label", (std::string) "lin_ow");
        amrex::ParmParse ppow("OceanWaves.lin_ow");
        ppow.add("type", (std::string) "LinearWaves");
        ppow.add("wave_height", 0.05);
        ppow.add("wave_length", 1.0);
        ppow.add("water_depth", 1.0);
        // Wave generation and numerical beach
        ppow.add("relax_zone_gen_length", 2.0);
        ppow.add("numerical_beach_length", 4.0);
    }
    {
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", 0.1);
    }
    {
        // Boundary conditions
        amrex::ParmParse ppxlo("xlo");
        ppxlo.add("type", (std::string) "wave_generation");
        ppxlo.addarr("velocity", gas_vel);
        ppxlo.add("vof", 0.0);
        ppxlo.add("density", 1.0);
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
    auto& repo = sim().repo();
    auto& velocity = repo.get_field("velocity");
    velocity.setVal(gas_vel, 3);
    oceanwaves.pre_init_actions();
    for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
        oceanwaves.initialize_fields(lev, mesh().Geom(lev));
    }
    // Do post-init step, which includes fillpatch calls
    oceanwaves.post_init_actions();

    auto& multiphase = sim().physics_manager().get<amr_wind::MultiPhase>();
    const amrex::Real rho1 = multiphase.rho1();
    const amrex::Real rho2 = multiphase.rho2();

    // Get fields for comparison
    auto& ow_vof = repo.get_field("ow_vof");
    auto& vof = repo.get_field("vof");
    auto& ow_velocity = repo.get_field("ow_velocity");
    auto& density = repo.get_field("density");

    // Check velocity field to confirm not modified
    amrex::Real error_total = bdy_error(vof, ow_vof);
    EXPECT_NEAR(error_total, 0.0, tol);
    make_target_velocity(ow_velocity, velocity, ow_vof);
    error_total = bdy_error(velocity, ow_velocity);
    EXPECT_NEAR(error_total, 0.0, tol);
    make_target_density(ow_vof, rho1, rho2);
    error_total = bdy_error(density, ow_vof);
    EXPECT_NEAR(error_total, 0.0, tol);
}

} // namespace amr_wind_tests
