#include "gtest/gtest.h"
#include "aw_test_utils/MeshTest.H"
#include "amr-wind/turbulence/TurbulenceModel.H"
#include "aw_test_utils/test_utils.H"

namespace amr_wind_tests {

namespace {

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
                const amrex::Real x = problo[0] + (i + offset) * dx[0];
                const amrex::Real y = problo[1] + (j + offset) * dx[1];
                const amrex::Real z = problo[2] + (k + offset) * dx[2];

                farr(i, j, k, 0) = x / sqrt(6.0) * srate;
                farr(i, j, k, 1) = y / sqrt(6.0) * srate;
                farr(i, j, k, 2) = z / sqrt(6.0) * srate;
            });
        }
    }
}

void init_field_amd(amr_wind::Field& fld, amrex::Real scale)
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
                const amrex::Real x = problo[0] + (i + offset) * dx[0];
                const amrex::Real y = problo[1] + (j + offset) * dx[1];
                const amrex::Real z = problo[2] + (k + offset) * dx[2];

                farr(i, j, k, 0) = 1 * x / sqrt(6.0) * scale;
                farr(i, j, k, 1) = -2 * y / sqrt(6.0) * scale;
                farr(i, j, k, 2) = -1 * z / sqrt(6.0) * scale;
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

class TurbLESTest : public MeshTest
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
        }
    }

    const amrex::Real dx = 10.0 / 10.0;
    const amrex::Real dy = 10.0 / 20.0;
    const amrex::Real dz = 10.0 / 30.0;
};

TEST_F(TurbLESTest, test_smag_setup_calc)
{
    // Parser inputs for turbulence model
    const amrex::Real Cs = 0.16;
    const amrex::Real visc = 1e-5;
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", (std::string) "Smagorinsky");
    }
    {
        amrex::ParmParse pp("Smagorinsky_coeffs");
        pp.add("Cs", Cs);
    }
    {
        amrex::ParmParse pp("transport");
        pp.add("viscosity", visc);
    }

    // Initialize necessary parts of solver
    populate_parameters();
    initialize_mesh();
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();

    // Create turbulence model
    sim().create_turbulence_model();
    // Get turbulence model
    auto& tmodel = sim().turbulence_model();

    // Get coefficients
    auto model_dict = tmodel.model_coeffs();

    for (const std::pair<const std::string, const amrex::Real> n : model_dict) {
        // Only a single model parameter, Cs
        EXPECT_EQ(n.first, "Cs");
        EXPECT_EQ(n.second, Cs);
    }

    // Constants for fields
    const amrex::Real srate = 0.5;
    const amrex::Real rho0 = 1.2;

    // Set up velocity field with constant strainrate
    auto& vel = sim().repo().get_field("velocity");
    init_field3(vel, srate);
    // Set up uniform unity density field
    auto& dens = sim().repo().get_field("density");
    dens.setVal(rho0);

    // Update turbulent viscosity directly
    tmodel.update_turbulent_viscosity(amr_wind::FieldState::New);
    auto& muturb = sim().repo().get_field("mu_turb");

    // Check values of turbulent viscosity
    auto min_val = utils::field_min(muturb);
    auto max_val = utils::field_max(muturb);
    const amrex::Real tol = 1e-12;
    const amrex::Real smag_answer =
        rho0 * std::pow(Cs, 2) * std::pow(std::cbrt(dx * dy * dz), 2) * srate;
    EXPECT_NEAR(min_val, smag_answer, tol);
    EXPECT_NEAR(max_val, smag_answer, tol);

    // Check values of effective viscosity
    auto& mueff = sim().repo().get_field("velocity_mueff");
    tmodel.update_mueff(mueff);
    min_val = utils::field_min(mueff);
    max_val = utils::field_max(mueff);
    EXPECT_NEAR(min_val, smag_answer + 1e-5, tol);
    EXPECT_NEAR(max_val, smag_answer + 1e-5, tol);

    // Check that this effective viscosity is what gets to icns diffusion
    auto visc_name = pde_mgr.icns().fields().mueff.name();
    EXPECT_EQ(visc_name, "velocity_mueff");
}

TEST_F(TurbLESTest, test_1eqKsgs_setup_calc)
{
    // Parser inputs for turbulence model
    const amrex::Real Ceps = 0.11;
    const amrex::Real Ce = 0.99;
    const amrex::Real Tref = 263.5;
    const amrex::Real gravz = 10.0;
    const amrex::Real rho0 = 1.2;
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
        pp.add("surface_temp_rate", -0.25);
        amrex::Vector<amrex::Real> t_hts{0.0, 100.0, 400.0};
        pp.addarr("temperature_heights", t_hts);
        amrex::Vector<amrex::Real> t_vals{265.0, 265.0, 268.0};
        pp.addarr("temperature_values", t_vals);
    }

    // Initialize necessary parts of solver
    populate_parameters();
    initialize_mesh();
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();

    // Create turbulence model
    sim().create_turbulence_model();
    // Get turbulence model
    auto& tmodel = sim().turbulence_model();

    // Get coefficients
    auto model_dict = tmodel.model_coeffs();

    int ct = 0;
    for (const std::pair<const std::string, const amrex::Real> n : model_dict) {
        // Two model parameters
        if (ct == 0) {
            EXPECT_EQ(n.first, "Ceps");
            EXPECT_EQ(n.second, Ceps);
        } else {
            EXPECT_EQ(n.first, "Ce");
            EXPECT_EQ(n.second, Ce);
        }
        ++ct;
    }

    // Constants for fields
    const amrex::Real srate = 20.0;
    const amrex::Real Tgz = 2.0;
    const amrex::Real tlscale_val = 1.1;
    const amrex::Real tke_val = 0.1;
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

    // Update turbulent viscosity directly
    tmodel.update_turbulent_viscosity(amr_wind::FieldState::New);
    auto& muturb = sim().repo().get_field("mu_turb");

    // Check values of turbulent viscosity
    const auto min_val = utils::field_min(muturb);
    const auto max_val = utils::field_max(muturb);
    const amrex::Real tol = 1e-12;
    const amrex::Real ksgs_answer =
        rho0 * Ce *
        amrex::min<amrex::Real>(
            std::cbrt(dx * dy * dz),
            0.76 * sqrt(tke_val / (Tgz * gravz) * Tref)) *
        sqrt(tke_val);
    EXPECT_NEAR(min_val, ksgs_answer, tol);
    EXPECT_NEAR(max_val, ksgs_answer, tol);
}

TEST_F(TurbLESTest, test_AMD_setup_calc)
{
    // Parser inputs for turbulence model
    const amrex::Real C = 0.3;
    const amrex::Real Tref = 200;
    const amrex::Real gravz = 10.0;
    const amrex::Real rho0 = 1.0;
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", (std::string) "AMD");
    }
    {
        amrex::ParmParse pp("AMD_coeffs");
        pp.add("C_poincare", C);
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
        amrex::Vector<amrex::Real> t_vals{200.0, 200.0, 200.0};
        pp.addarr("temperature_values", t_vals);
    }

    // Initialize necessary parts of solver
    populate_parameters();
    initialize_mesh();
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();

    // Create turbulence model
    sim().create_turbulence_model();
    // Get turbulence model
    auto& tmodel = sim().turbulence_model();

    // Get coefficients
    auto model_dict = tmodel.model_coeffs();

    for (const std::pair<const std::string, const amrex::Real> n : model_dict) {
        // Only a single model parameter, Cs
        EXPECT_EQ(n.first, "C_poincare");
        EXPECT_EQ(n.second, C);
    }

    // Constants for fields
    const amrex::Real scale = 1.50;
    const amrex::Real Tgz = 20.0;
    // Set up velocity field with constant strainrate
    auto& vel = sim().repo().get_field("velocity");
    init_field_amd(vel, scale);
    // Set up uniform unity density field
    auto& dens = sim().repo().get_field("density");
    dens.setVal(rho0);
    // Set up temperature field with constant gradient in z
    auto& temp = sim().repo().get_field("temperature");
    init_field1(temp, Tgz);

    // Update turbulent viscosity directly
    tmodel.update_turbulent_viscosity(amr_wind::FieldState::New);
    auto& muturb = sim().repo().get_field("mu_turb");

    // Check values of turbulent viscosity
    const auto min_val = utils::field_min(muturb);
    const auto max_val = utils::field_max(muturb);
    const amrex::Real tol = 1e-12;

    const amrex::Real amd_answer =
        C *
        (-2.0 * std::pow(scale / sqrt(6), 3) *
             (dx * dx - 8 * dy * dy - dz * dz) +
         gravz / Tref * (-1.0 * Tgz * scale / sqrt(6) * dz * dz)) /
        (1 * scale * scale);
    EXPECT_NEAR(min_val, amd_answer, tol);
    EXPECT_NEAR(max_val, amd_answer, tol);

    // Check values of alphaeff
    auto& alphaeff = sim().repo().declare_cc_field("alphaeff");
    tmodel.update_alphaeff(alphaeff);
    const auto ae_min_val = utils::field_min(alphaeff);
    const auto ae_max_val = utils::field_max(alphaeff);
    const amrex::Real amd_ae_answer = C * dz * dz * scale * 1.0 / sqrt(6);
    EXPECT_NEAR(ae_min_val, amd_ae_answer, tol);
    EXPECT_NEAR(ae_max_val, amd_ae_answer, tol);
}
} // namespace amr_wind_tests
