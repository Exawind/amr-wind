#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"

#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/equation_systems/icns/icns.H"
#include "amr-wind/equation_systems/icns/icns_ops.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/equation_systems/icns/source_terms/BodyForce.H"
#include "amr-wind/equation_systems/icns/source_terms/ABLForcing.H"
#include "amr-wind/equation_systems/icns/source_terms/GeostrophicForcing.H"
#include "amr-wind/equation_systems/icns/source_terms/CoriolisForcing.H"

namespace {
void write_target_velocity_file(const std::string& fname)
{
    std::ofstream os(fname);
    // Write header **??**
    os << "Time\tWindSpeed\tHorzAngle\n";
    // Write time table
    os << "0.0\t8.0\t0.0\n";
    os << "0.2\t8.0\t5.0\n";
    os << "5.0\t8.0\t30.0\n";
}
} // namespace

namespace amr_wind_tests {

using ICNSFields =
    amr_wind::pde::FieldRegOp<amr_wind::pde::ICNS, amr_wind::fvm::Godunov>;

// Tests in this file involve ABLForcing, BodyForce, and GeostrophicWind
class ABLSrcTimeTableTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        // Make computational domain like ABL mesh
        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{8, 8, 64}};
            pp.addarr("n_cell", ncell);
        }

        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> probhi{{120.0, 120.0, 1000.0}};
            pp.addarr("prob_hi", probhi);
        }

        // Parameters for initializing ABL case
        {
            amrex::ParmParse pp("ABL");
            amrex::Vector<amrex::Real> theights{{0.0, 650.0, 750.0, 1000.0}};
            amrex::Vector<amrex::Real> tvalues{{300.0, 300.0, 308.0, 308.75}};
            pp.addarr("temperature_heights", theights);
            pp.addarr("temperature_values", tvalues);
            pp.add("perturb_ref_height", 50.0);
            pp.add("reference_temperature", 300.0);
            pp.add("kappa", 0.41);
            pp.add("surface_roughness_z0", 0.1);
        }

        {
            // Physics
            amrex::ParmParse pp("incflo");
            amrex::Vector<std::string> phystr{"ABL"};
            pp.addarr("physics", phystr);
            pp.add("density", 1.0);
        }

        // Timestep size
        {
            amrex::ParmParse pp("time");
            pp.add("fixed_dt", 0.1);
        }
    }
    std::string tvel_fname = "target_velocities.txt";
    std::string forces_fname = "abl_forces.txt";
    amrex::Real dt{0.1};
};

TEST_F(ABLSrcTimeTableTest, abl)
{
    constexpr amrex::Real tol = 1.0e-12;

    // Write target wind file
    write_target_velocity_file(tvel_fname);

    // Set up simulation parameters and mesh
    populate_parameters();
    // ABL Forcing
    {
        amrex::ParmParse pp("ABLForcing");
        pp.add("abl_forcing_height", 90.0);
        pp.add("velocity_timetable", tvel_fname);
        pp.add("forcing_timetable_output_file", forces_fname);
    }
    initialize_mesh();

    // Set up PDEs and physics objects
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();
    auto& velocity = pde_mgr.icns().fields().field;
    auto& ABL = sim().physics_manager().get<amr_wind::ABL>();

    // Get icns source term, which will be tested
    auto& src_term = pde_mgr.icns().fields().src_term;
    src_term.setVal(0.0);
    // Source term object for abl_forcing
    amr_wind::pde::icns::ABLForcing abl_forcing(sim());
    run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.tilebox();
        const auto& src_arr = src_term(lev).array(mfi);

        abl_forcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
    });

    // Initial velocity should be the same as target velocity, so force is 0
    amrex::Vector<amrex::Real> target_force = {0.0, 0.0, 0.0};
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        const auto min_val = utils::field_min(src_term, i);
        const auto max_val = utils::field_max(src_term, i);
        EXPECT_NEAR(min_val, target_force[i], tol);
        EXPECT_NEAR(min_val, max_val, tol);
    }

    // Change mean velocity and recalculate
    const amrex::Vector<amrex::Real> init_vel{8.0, 0.0, 0.0};
    const amrex::Vector<amrex::Real> new_vel{5.0, 3.0, 0.0};
    // This function actually calculates the forcing terms, applied later
    // This function also is the one that writes forces to file
    abl_forcing.set_mean_velocities(new_vel[0], new_vel[1]);
    // Doing this through ABL physics is complicated and would require
    // calc_averages in ABLStats, then pre_advance_work in ABL physics
    src_term.setVal(0.0);
    run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.tilebox();
        const auto& src_arr = src_term(lev).array(mfi);

        abl_forcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
    });
    // Velocity at hub height is (5, 3) and target is (8, 0)
    target_force[0] = (init_vel[0] - new_vel[0]) / dt;
    target_force[1] = (init_vel[1] - new_vel[1]) / dt;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        const auto min_val = utils::field_min(src_term, i);
        const auto max_val = utils::field_max(src_term, i);
        EXPECT_NEAR(min_val, target_force[i], tol);
        EXPECT_NEAR(min_val, max_val, tol);
    }

    // Advance time (twice to make current_time change)
    // sim().time().set_current_cfl(0.1, 0.1, 0.1);
    sim().time().new_timestep();
    sim().time().new_timestep();

    // Go back to original mean values on mesh
    abl_forcing.set_mean_velocities(init_vel[0], init_vel[1]);
    // Recalculate forcing and check
    src_term.setVal(0.0);
    run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.tilebox();
        const auto& src_arr = src_term(lev).array(mfi);

        abl_forcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
    });
    // Velocity at hub height is 8 at 0deg and target is 8 at 2.5deg
    target_force[0] = (8.0 * std::cos(M_PI / 180.0 * 2.5) - init_vel[0]) / dt;
    target_force[1] = (8.0 * std::sin(M_PI / 180.0 * 2.5) - init_vel[1]) / dt;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        const auto min_val = utils::field_min(src_term, i);
        const auto max_val = utils::field_max(src_term, i);
        EXPECT_NEAR(min_val, target_force[i], tol);
        EXPECT_NEAR(min_val, max_val, tol);
    }

    // Delete target wind file
    const char* fname = tvel_fname.c_str();
    {
        std::ifstream f(fname);
        if (f.good()) {
            remove(fname);
        }
        // Check that file is removed
        std::ifstream ff(fname);
        EXPECT_FALSE(ff.good());
    }
}

TEST_F(ABLSrcTimeTableTest, bodyforce)
{

    // Body Forcing
    {
        amrex::ParmParse pp("BodyForce");
        pp.add("type", "uniform_timetable");
        pp.add("uniform_timetable_file", forces_fname);
    }
}

TEST_F(ABLSrcTimeTableTest, geostrophic)
{
    constexpr amrex::Real tol = 1.0e-12;

    // write target wind file

    // Geostrophic Forcing
    {
        amrex::ParmParse pp("GeostrophicForcing");
        pp.add("geostrophic_wind_timetable", tvel_fname);
    }

    // initialize and check force values

    // check forces at later time

    // delete target wind file
}

} // namespace amr_wind_tests
