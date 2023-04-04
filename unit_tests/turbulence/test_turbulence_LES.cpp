#include "gtest/gtest.h"
#include "aw_test_utils/MeshTest.H"
#include "amr-wind/turbulence/TurbulenceModel.H"
#include "aw_test_utils/test_utils.H"

namespace amr_wind_tests {

namespace {

void init_field(amr_wind::Field& fld, amrex::Real srate)
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

TEST_F(TurbLESTest, test_turb_smag_setup)
{
    // Parser inputs for turbulence model
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", (std::string) "Smagorinsky");
    }
    {
        amrex::ParmParse pp("Smagorinsky_coeffs");
        pp.add("Cs", 0.16);
    }
    {
        amrex::ParmParse pp("transport");
        pp.add("viscosity", 1e-5);
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
        EXPECT_EQ(n.second, 0.16);
    }

    // Set up velocity field with constant strainrate
    auto& vel = sim().repo().get_field("velocity");
    init_field(vel, 1.0);
    // Set up uniform unity density field
    auto& dens = sim().repo().get_field("density");
    dens.setVal(1.0);

    // Update turbulent viscosity by updating effective viscosity
    // auto& mueff = sim().repo().get_field("velocity_mueff");
    // tmodel.update_mueff(mueff);
    tmodel.update_turbulent_viscosity(amr_wind::FieldState::New);
    auto& muturb = sim().repo().get_field("mu_turb");

    // Check values of turbulent viscosity
    auto min_val = utils::field_min(muturb);
    auto max_val = utils::field_max(muturb);
    const amrex::Real tol = 1e-12;
    const amrex::Real smag_answer =
        1.0 * std::pow(0.16, 2) * std::pow(std::cbrt(dx * dy * dz), 2);
    EXPECT_NEAR(min_val, smag_answer, tol);
    EXPECT_NEAR(max_val, smag_answer, tol);

    // Check values of effective viscosity
    auto& mueff = sim().repo().get_field("velocity_mueff");
    tmodel.update_mueff(mueff);
    min_val = utils::field_min(mueff);
    max_val = utils::field_max(mueff);
    EXPECT_NEAR(min_val, smag_answer + 1e-5, tol);
    EXPECT_NEAR(max_val, smag_answer + 1e-5, tol);
}

} // namespace amr_wind_tests
