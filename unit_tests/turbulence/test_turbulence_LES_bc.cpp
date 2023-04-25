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
                const amrex::Real x = problo[0] + (i + offset) * dx[0];
                const amrex::Real y = problo[1] + (j + offset) * dx[1];
                const amrex::Real z = problo[2] + (k + offset) * dx[2];

                farr(i, j, k, 0) = z / sqrt(2.0) * srate;
                farr(i, j, k, 1) = z / sqrt(2.0) * srate;
                farr(i, j, k, 2) = 0.0;
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
        }
        {
            // Periodicity
            amrex::ParmParse pp("geometry");
            amrex::Vector<int> periodic{{1, 1, 0}};
            pp.addarr("is_periodic", periodic);
            // Boundary conditions
            amrex::ParmParse ppzhi("zhi");
            ppzhi.add("type", (std::string) "slip_wall");
            // zlo is defined in each case
        }
    }

    const amrex::Real dx = 10.0 / 10.0;
    const amrex::Real dy = 10.0 / 20.0;
    const amrex::Real dz = 10.0 / 30.0;
    const amrex::Real tol = 1.0e-12;
};

TEST_F(TurbLESTestBC, test_1eqKsgs_noslip)
{
    // Parser inputs for turbulence model
    const amrex::Real Ceps = 0.11;
    const amrex::Real Ce = 0.99;
    const amrex::Real Tref = 263.5;
    const amrex::Real gravz = 10.0;
    const amrex::Real rho0 = 1.2;
    {
        amrex::ParmParse pp("zlo");
        pp.add("type", (std::string) "no_slip_wall");
    }
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

    // Constants for fields
    const amrex::Real srate = 20.0;
    const amrex::Real Tgz = 2.0;
    const amrex::Real tlscale_val = 1.1;
    const amrex::Real tke_val = 0.1;
    // Set up velocity field with constant strainrate
    auto& vel = sim().repo().get_field("velocity");
    init_field3(vel, srate);
    // Perform fillpatch to allow BCs to operate
    vel.fillpatch(0.0);
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

    // Get shear production field
    auto& shear_prod = sim().repo().get_field("shear_prod");

    // Check for constant value in bulk
    auto shear_bulk = get_val_at_kindex(shear_prod, muturb, 0, 1) / 10. / 20.;
    EXPECT_NEAR(shear_bulk, srate, tol);

    // Velocity gradients assumed in setup (init_field3)
    const amrex::Real uz_bulk = srate / sqrt(2.0);
    const amrex::Real uz_wallcell = (uz_bulk * 1.5 * dz - 0.0) / (2.0 * dz);
    const amrex::Real uz_wallface =
        (0.5 * (uz_bulk * 1.5 * dz + uz_bulk * 0.5 * dz) - 0.0) / dz;
    // Naive cell-centered answer, with no_slip_wall (Dirichlet)
    const amrex::Real s_naive = uz_wallcell * sqrt(2.0);
    // Check for different value just above wall due to BC
    auto shear_wall = get_val_at_kindex(shear_prod, muturb, 0, 0) / 10. / 20.;
    // Check that the result is not equal to the naive value
    EXPECT_GT(std::abs(shear_wall - s_naive), tol);
    // Answer that accounts for location of wall at cell face
    const amrex::Real s_true = uz_wallface * sqrt(2.0);
    EXPECT_NEAR(shear_wall, s_true, tol);
}

} // namespace amr_wind_tests
