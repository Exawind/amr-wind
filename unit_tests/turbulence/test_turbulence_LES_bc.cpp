#include "gtest/gtest.h"
#include "aw_test_utils/MeshTest.H"
#include "amr-wind/turbulence/TurbulenceModel.H"
#include "aw_test_utils/test_utils.H"

namespace amr_wind_tests {

namespace {

amrex::Real get_val_at_kindex(
    amr_wind::Field& field,
    amr_wind::Field& divisor,
    const int comp,
    const int kref)
{
    const int lev = 0;
    amrex::Real error_total = 0;

    error_total += amrex::ReduceSum(
        field(lev), divisor(lev), 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx, amrex::Array4<amrex::Real const> const& f_arr,
            amrex::Array4<amrex::Real const> const& div_arr) -> amrex::Real {
            amrex::Real error = 0;

            amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                // Check if current cell is just above lower wall
                if (k == kref) {
                    // Add field value to output
                    error += sqrt(f_arr(i, j, k, comp) / div_arr(i, j, k));
                }
            });

            return error;
        });
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    return error_total;
}

void init_field3(amr_wind::Field& fld, amrex::Real srate)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();

    amrex::Real offset = 0.0;
    if (fld.field_location() == amr_wind::FieldLoc::CELL) {
        offset = 0.5;
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();

        for (amrex::MFIter mfi(fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.growntilebox();
            const auto& farr = fld(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const amrex::Real z = problo[2] + (k + offset) * dx[2];

                farr(i, j, k, 0) = z / 2.0 * srate + 2.0;
                farr(i, j, k, 1) = z / 2.0 * srate + 2.0;
                farr(i, j, k, 2) = z / 2.0 * srate + 2.0;
            });
        }
    }
}

void init_field1(amr_wind::Field& fld, amrex::Real tgrad)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();

    amrex::Real offset = 0.0;
    if (fld.field_location() == amr_wind::FieldLoc::CELL) {
        offset = 0.5;
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();

        for (amrex::MFIter mfi(fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.growntilebox();
            const auto& farr = fld(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const amrex::Real z = problo[2] + (k + offset) * dx[2];

                farr(i, j, k, 0) = z * tgrad;
            });
        }
    }
}

} // namespace

class TurbLESTestBC : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{10, 20, 30}};
            pp.addarr("n_cell", ncell);
            pp.add("blocking_factor", 2);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
            amrex::Vector<amrex::Real> probhi{{10.0, 10.0, 10.0}};
            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
            amrex::Vector<int> periodic{{1, 1, 0}};
            pp.addarr("is_periodic", periodic);
            // Boundary conditions
            amrex::ParmParse ppzhi("zhi");
            ppzhi.add("type", (std::string) "slip_wall");
            // zlo is defined in each case
        }
    }
    void OneEqKsgs_setup_params() const
    {
        {
            amrex::ParmParse pp("turbulence");
            pp.add("model", (std::string) "OneEqKsgsM84");
        }
        {
            amrex::ParmParse pp("OneEqKsgsM84_coeffs");
            pp.add("Ceps", Ceps);
            pp.add("Ce", Ce);
        }
        {
            amrex::ParmParse pp("incflo");
            amrex::Vector<std::string> physics{"ABL"};
            pp.addarr("physics", physics);
            pp.add("density", rho0);
            amrex::Vector<amrex::Real> vvec{8.0, 0.0, 0.0};
            pp.addarr("velocity", vvec);
            amrex::Vector<amrex::Real> gvec{0.0, 0.0, -gravz};
            pp.addarr("gravity", gvec);
        }
        {
            amrex::ParmParse pp("ABL");
            pp.add("reference_temperature", Tref);
            amrex::Vector<amrex::Real> t_hts{0.0, 100.0, 400.0};
            pp.addarr("temperature_heights", t_hts);
            amrex::Vector<amrex::Real> t_vals{265.0, 265.0, 268.0};
            pp.addarr("temperature_values", t_vals);
        }
        {
            amrex::ParmParse pp("transport");
            pp.add("viscosity", mu);
        }
    }
    void test_calls_body(bool do_postsolve = false)
    {
        // Initialize necessary parts of solver
        initialize_mesh();
        auto& pde_mgr = sim().pde_manager();
        auto& icns_eq = pde_mgr.register_icns();
        sim().init_physics();
        sim().create_turbulence_model();
        icns_eq.initialize();

        // Get turbulence model
        auto& tmodel = sim().turbulence_model();
        // Set up velocity field with constant strainrate
        auto& vel = sim().repo().get_field("velocity");
        init_field3(vel, srate);
        // Set up uniform unity density field
        auto& dens = sim().repo().get_field("density");
        dens.setVal(rho0);
        // Set up temperature field with constant gradient in z
        auto& temp = sim().repo().get_field("temperature");
        init_field1(temp, Tgz);
        // Give values to tlscale and tke arrays
        auto& tlscale = sim().repo().get_field("turb_lscale");
        tlscale.setVal(tlscale_val);
        auto& tke = sim().repo().get_field("tke");
        tke.setVal(tke_val);

        // Perform post init for physics
        // (turns on wall model if specified)
        for (auto& pp : sim().physics()) {
            pp->post_init_actions();
        }
        // Perform fillpatch to allow BCs to operate
        vel.fillpatch(0.0);
        temp.fillpatch(0.0);
        // Ensure nonzero viscosity (only molecular at this point)
        icns_eq.compute_mueff(amr_wind::FieldState::New);

        // Advance states (important for wall model)
        pde_mgr.advance_states();
        // Compute diffusion term to use BCs (wall model)
        icns_eq.compute_diffusion_term(amr_wind::FieldState::New);

        // Post_solve afterward
        if (do_postsolve) {
            icns_eq.post_solve_actions();
        }

        // Update turbulent viscosity directly
        tmodel.update_turbulent_viscosity(amr_wind::FieldState::New);
    }

    const amrex::Real dx = 10.0 / 10.0;
    const amrex::Real dy = 10.0 / 20.0;
    const amrex::Real dz = 10.0 / 30.0;
    const amrex::Real tol = 1.0e-12;
    // Parser inputs for turbulence model
    const amrex::Real Ceps = 0.11;
    const amrex::Real Ce = 0.99;
    const amrex::Real Tref = 263.5;
    const amrex::Real gravz = 10.0;
    const amrex::Real rho0 = 1.2;
    // Constants for fields
    const amrex::Real srate = 20.0;
    const amrex::Real Tgz = 2.0;
    const amrex::Real tlscale_val = 1.1;
    const amrex::Real tke_val = 0.1;
    // Molecular viscosity
    const amrex::Real mu = 1e-2;
};

TEST_F(TurbLESTestBC, test_1eqKsgs_noslip)
{
    populate_parameters();
    {
        amrex::ParmParse pp("zlo");
        pp.add("type", (std::string) "no_slip_wall");
    }
    OneEqKsgs_setup_params();
    test_calls_body();
    auto& muturb = sim().repo().get_field("mu_turb");

    // Get shear production field
    auto& shear_prod = sim().repo().get_field("shear_prod");

    // Check for constant value in bulk
    auto shear_bulk = get_val_at_kindex(shear_prod, muturb, 0, 1) / 10. / 20.;
    EXPECT_NEAR(shear_bulk, srate, tol);

    // Velocity gradients assumed in setup (init_field3)
    const amrex::Real uz_bulk = srate / 2.0;

    // Naive cell-centered answer, with no_slip_wall (Dirichlet)
    const amrex::Real uz_wallcell =
        (uz_bulk * 1.5 * dz + 2.0 - 0.0) / (2.0 * dz);
    // Wall-normal direction
    const amrex::Real wz_wallcell = uz_wallcell;

    const amrex::Real s_naive =
        sqrt(2.0 * wz_wallcell * wz_wallcell + 2.0 * uz_wallcell * uz_wallcell);
    // Check for different value just above wall due to BC
    auto shear_wall = get_val_at_kindex(shear_prod, muturb, 0, 0) / 10. / 20.;
    // Check that the result is not equal to the naive value
    EXPECT_GT(std::abs(shear_wall - s_naive), tol);
    // Answer that accounts for location of wall at cell face
    const amrex::Real uz_wallface =
        ((uz_bulk * 1.5 * dz + 2.0) / 3.0 + uz_bulk * 0.5 * dz + 2.0 - 0.0) /
        dz;
    const amrex::Real wz_wallface = uz_wallface;
    const amrex::Real s_true =
        sqrt(2.0 * wz_wallface * wz_wallface + 2.0 * uz_wallface * uz_wallface);
    EXPECT_NEAR(shear_wall, s_true, tol);
}

TEST_F(TurbLESTestBC, test_1eqKsgs_slip)
{
    populate_parameters();
    {
        amrex::ParmParse pp("zlo");
        pp.add("type", (std::string) "slip_wall");
    }
    OneEqKsgs_setup_params();
    test_calls_body();
    auto& muturb = sim().repo().get_field("mu_turb");

    // Get shear production field
    auto& shear_prod = sim().repo().get_field("shear_prod");

    // Check for constant value in bulk
    auto shear_bulk = get_val_at_kindex(shear_prod, muturb, 0, 1) / 10. / 20.;
    EXPECT_NEAR(shear_bulk, srate, tol);

    // Check for different value just above wall due to BC
    auto shear_wall = get_val_at_kindex(shear_prod, muturb, 0, 0) / 10. / 20.;
    // Answer that accounts for slip_wall BC
    // (tangential extrapolation, dirichlet normal)
    const amrex::Real uz_bulk = srate / 2.0;
    const amrex::Real wz_wall =
        ((uz_bulk * 1.5 * dz + 2.0) / 3.0 + uz_bulk * 0.5 * dz + 2.0 - 0.0) /
        dz;
    const amrex::Real s_true =
        sqrt(2.0 * wz_wall * wz_wall + 2.0 * uz_bulk * uz_bulk);
    EXPECT_NEAR(shear_wall, s_true, tol);
}

TEST_F(TurbLESTestBC, test_1eqKsgs_wallmodel)
{
    constexpr amrex::Real kappa = 0.4;
    constexpr amrex::Real z0 = 0.11;
    populate_parameters();
    {
        amrex::ParmParse pp("zlo");
        pp.add("type", (std::string) "wall_model");
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("wall_shear_stress_type", (std::string) "local");
        pp.add("kappa", kappa);
        pp.add("surface_roughness_z0", z0);
    }
    {
        // Explicit diffusion populates wall cells, doesn't update velocity
        amrex::ParmParse pp("incflo");
        pp.add("diffusion_type", 0);
    }
    OneEqKsgs_setup_params();
    const bool do_postsolve = true;
    test_calls_body(do_postsolve);
    auto& muturb = sim().repo().get_field("mu_turb");

    // Get shear production field
    auto& shear_prod = sim().repo().get_field("shear_prod");

    // Check for constant value in bulk
    auto shear_bulk = get_val_at_kindex(shear_prod, muturb, 0, 1) / 10. / 20.;
    EXPECT_NEAR(shear_bulk, srate, tol);

    // Velocity gradients assumed in setup (init_field3)
    const amrex::Real uz_bulk = srate / 2.0;
    // Get value just above wall due to BC
    auto shear_wall = get_val_at_kindex(shear_prod, muturb, 0, 0) / 10. / 20.;
    // (tangential extrapolation, dirichlet normal)
    const amrex::Real wz_wall =
        ((uz_bulk * 1.5 * dz + 2.0) / 3.0 + uz_bulk * 0.5 * dz + 2.0 - 0.0) /
        dz;
    const amrex::Real s_true =
        sqrt(2.0 * wz_wall * wz_wall + 2.0 * uz_bulk * uz_bulk);
    EXPECT_NEAR(shear_wall, s_true, tol);
}

TEST_F(TurbLESTestBC, test_1eqKsgs_wallmodel_failnofillpatch)
{
    constexpr amrex::Real kappa = 0.4;
    constexpr amrex::Real z0 = 0.11;
    populate_parameters();
    {
        amrex::ParmParse pp("zlo");
        pp.add("type", (std::string) "wall_model");
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("wall_shear_stress_type", (std::string) "local");
        pp.add("kappa", kappa);
        pp.add("surface_roughness_z0", z0);
    }
    {
        // Explicit diffusion populates wall cells, doesn't update velocity
        amrex::ParmParse pp("incflo");
        pp.add("diffusion_type", 0);
    }
    OneEqKsgs_setup_params();
    const bool do_postsolve = false;
    test_calls_body(do_postsolve);
    auto& muturb = sim().repo().get_field("mu_turb");

    // Get shear production field
    auto& shear_prod = sim().repo().get_field("shear_prod");

    // Check for constant value in bulk
    auto shear_bulk = get_val_at_kindex(shear_prod, muturb, 0, 1) / 10. / 20.;
    EXPECT_NEAR(shear_bulk, srate, tol);

    // Velocity gradients assumed in setup (init_field3)
    const amrex::Real uz_bulk = srate / 2.0;

    // Get value just above wall due to BC
    auto shear_wall = get_val_at_kindex(shear_prod, muturb, 0, 0) / 10. / 20.;
    // (tangential extrapolation, dirichlet normal)
    const amrex::Real wz_wall =
        ((uz_bulk * 1.5 * dz + 2.0) / 3.0 + uz_bulk * 0.5 * dz + 2.0 - 0.0) /
        dz;
    const amrex::Real s_true =
        sqrt(2.0 * wz_wall * wz_wall + 2.0 * uz_bulk * uz_bulk);
    // This is checking the correct value -- without the fillpatch, it is wrong
    // (passing test indicates fillpatch is needed after diffusion performed)
    EXPECT_GT(std::abs(shear_wall - s_true), tol);

    // Shear value at the wall used by wall model
    const amrex::Real zref = 0.5 * dz;
    const amrex::Real uref = zref * uz_bulk + 2.0;
    const amrex::Real vmag_ref = std::sqrt(2.0 * uref * uref);
    const amrex::Real utau = kappa * vmag_ref / (std::log(zref / z0));
    const amrex::Real uz_wm = uref / vmag_ref * std::pow(utau, 2) * rho0 / mu;

    // Velocity gradient with wallmodel value included as dirichlet
    const amrex::Real uz_wmdirichlet =
        (1.0 / 3.0 * (uz_bulk * 1.5 * dz + 2.0) +
         1.0 * (uz_bulk * 0.5 * dz + 2.0) - 4.0 / 3.0 * uz_wm) /
        dz;

    // Naive answer, assumes wall Dirichlet
    const amrex::Real s_naive =
        sqrt(2.0 * wz_wall * wz_wall + 2.0 * uz_wmdirichlet * uz_wmdirichlet);
    // This is checking the wrong value
    // (passing test indicates it is wrong in expected way w/out fillpatch)
    EXPECT_NEAR(shear_wall, s_naive, std::abs(shear_wall) * tol);
}

TEST_F(TurbLESTestBC, test_1eqKsgs_zerogradient)
{
    populate_parameters();
    {
        amrex::ParmParse pp("zlo");
        pp.add("type", (std::string) "zero_gradient");
    }
    OneEqKsgs_setup_params();
    test_calls_body();
    auto& muturb = sim().repo().get_field("mu_turb");

    // Get shear production field
    auto& shear_prod = sim().repo().get_field("shear_prod");

    // Check for constant value in bulk
    auto shear_bulk = get_val_at_kindex(shear_prod, muturb, 0, 1) / 10. / 20.;
    EXPECT_NEAR(shear_bulk, srate, tol);

    // Velocity gradients assumed in setup (init_field3)
    const amrex::Real uz_bulk = srate / 2.0;

    // Naive answer, assumes extrapolation
    const amrex::Real s_naive = srate;
    // Check for different value just above wall due to BC
    auto shear_wall = get_val_at_kindex(shear_prod, muturb, 0, 0) / 10. / 20.;
    // Check that the result is not equal to the naive value
    EXPECT_GT(std::abs(shear_wall - s_naive), tol);
    // Answer that accounts for BC (Neumann)
    const amrex::Real uz_wallface_neumann =
        (1.0 / 3.0 * (uz_bulk * 1.5 * dz + 2.0) +
         1.0 * (uz_bulk * 0.5 * dz + 2.0) -
         4.0 / 3.0 * (uz_bulk * 0.5 * dz + 2.0)) /
        dz;
    // Neumann is applied to w as well
    const amrex::Real wz_wallface_neumann = uz_wallface_neumann;
    const amrex::Real s_true = sqrt(
        2.0 * wz_wallface_neumann * wz_wallface_neumann +
        2.0 * uz_wallface_neumann * uz_wallface_neumann);
    EXPECT_NEAR(shear_wall, s_true, tol);
}

TEST_F(TurbLESTestBC, test_1eqKsgs_symmetricwall)
{
    populate_parameters();
    {
        amrex::ParmParse pp("zlo");
        pp.add("type", (std::string) "symmetric_wall");
    }
    OneEqKsgs_setup_params();
    test_calls_body();
    auto& muturb = sim().repo().get_field("mu_turb");

    // Get shear production field
    auto& shear_prod = sim().repo().get_field("shear_prod");

    // Check for constant value in bulk
    auto shear_bulk = get_val_at_kindex(shear_prod, muturb, 0, 1) / 10. / 20.;
    EXPECT_NEAR(shear_bulk, srate, tol);

    // Velocity gradients assumed in setup (init_field3)
    const amrex::Real uz_bulk = srate / 2.0;

    // Naive answer, assumes extrapolation
    const amrex::Real s_naive = srate;
    // Check for different value just above wall due to BC
    auto shear_wall = get_val_at_kindex(shear_prod, muturb, 0, 0) / 10. / 20.;
    // Check that the result is not equal to the naive value
    EXPECT_GT(std::abs(shear_wall - s_naive), tol);
    // Answer that accounts for BC (Neumann)
    const amrex::Real uz_wallface_neumann =
        (1.0 / 3.0 * (uz_bulk * 1.5 * dz + 2.0) +
         1.0 * (uz_bulk * 0.5 * dz + 2.0) -
         4.0 / 3.0 * (uz_bulk * 0.5 * dz + 2.0)) /
        dz;
    // Wall condition on w
    const amrex::Real wz_wallface =
        (1.0 / 3.0 * (uz_bulk * 1.5 * dz + 2.0) +
         1.0 * (uz_bulk * 0.5 * dz + 2.0) - 4.0 / 3.0 * (0.0)) /
        dz;
    const amrex::Real s_true = sqrt(
        2.0 * wz_wallface * wz_wallface +
        2.0 * uz_wallface_neumann * uz_wallface_neumann);
    EXPECT_NEAR(shear_wall, s_true, tol);
}

} // namespace amr_wind_tests
