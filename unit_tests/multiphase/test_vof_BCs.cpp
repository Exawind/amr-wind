#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/equation_systems/vof/vof.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"

namespace amr_wind_tests {

namespace {
void initialize_volume_fractions(
    const amrex::Box& bx, const amrex::Array4<amrex::Real>& vof_arr)
{
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        vof_arr(i, j, k) = 1.0 - 0.1 * (i + j + k);
    });
}
} // namespace

class VOFBCTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{2, 2, 2}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 2);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
            amrex::Vector<amrex::Real> probhi{{1.0, 1.0, 1.0}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
        {
            amrex::ParmParse pp("MultiPhase");
            pp.add("density_fluid1", m_rho1);
            pp.add("density_fluid2", m_rho2);
        }
        {
            amrex::ParmParse pp("incflo");
            amrex::Vector<std::string> physics{"MultiPhase"};
            pp.addarr("physics", physics);
            pp.add("use_godunov", (int)1);
            pp.add("godunov_type", (std::string) "weno");
        }
        {
            amrex::ParmParse pp("time");
            pp.add("fixed_dt", dt);
        }
    }

    void
    testing_BC_coorddir(const int option, const int dir, const amrex::Real tol)
    {

        // Get string version of direction
        std::string sdir = "d";
        std::string odir0 = "d";
        std::string odir1 = "d";
        switch (dir) {
        case 0:
            sdir = "x";
            odir0 = "y";
            odir1 = "z";
            break;
        case 1:
            sdir = "y";
            odir0 = "x";
            odir1 = "z";
            break;
        case 2:
            sdir = "z";
            odir0 = "x";
            odir1 = "y";
            break;
        }

        // Set up "correct" answers according to each option
        int vof_distr = 0;
        bool nonzero_flux = true;
        std::string bc_string = "sixteen chars --";
        switch (option) {
        case 1:
            // Mass inflow
            vof_distr = 0;       // prescribed vof value in bdy
            nonzero_flux = true; // flux can occur
            bc_string = "mass_inflow";
            break;
        case 2:
            // Slip wall
            vof_distr = 1;        // high-order extrapolation
            nonzero_flux = false; // flux cannot occur
            bc_string = "slip_wall";
            break;
        case 3:
            // No-slip wall
            vof_distr = 2;        // zero-gradient at bdy
            nonzero_flux = false; // flux cannot occur
            bc_string = "no_slip_wall";
            break;
        case 4:
            // Pressure outflow
            vof_distr = 2;       // zero-gradient at bdy
            nonzero_flux = true; // flux can occur
            bc_string = "pressure_outflow";
            break;
        }

        populate_parameters();
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<int> periodic{{0, 0, 0}};
            pp.addarr("is_periodic", periodic);
        }
        {
            amrex::ParmParse pp(sdir + "lo");
            pp.add("type", bc_string);
            if (option == 1) {
                // Specify vof
                pp.add("vof", m_vof_bdyval);
                // Specify other quantities to avoid errors, not actually used
                pp.add("levelset", 0.0);
                pp.add("density", 1.0);
                pp.addarr(
                    "velocity", amrex::Vector<amrex::Real>{0.0, 0.0, 0.0});
            }
        }
        {
            amrex::ParmParse pp(sdir + "hi");
            pp.add("type", bc_string);
            if (option == 1) {
                pp.add("vof", m_vof_bdyval);
                pp.add("levelset", 0.0);
                pp.add("density", 1.0);
                pp.addarr(
                    "velocity", amrex::Vector<amrex::Real>{0.0, 0.0, 0.0});
            }
        }
        {
            amrex::ParmParse pp(odir0 + "lo");
            pp.add("type", (std::string) "slip_wall");
        }
        {
            amrex::ParmParse pp(odir0 + "hi");
            pp.add("type", (std::string) "slip_wall");
        }
        {
            amrex::ParmParse pp(odir1 + "lo");
            pp.add("type", (std::string) "slip_wall");
        }
        {
            amrex::ParmParse pp(odir1 + "hi");
            pp.add("type", (std::string) "slip_wall");
        }

        initialize_mesh();

        auto& repo = sim().repo();

        // Set up icns PDE to set up mac velocity, access to VOF PDE
        auto& pde_mgr = sim().pde_manager();
        auto& mom_eqn = pde_mgr.register_icns();
        mom_eqn.initialize();

        // Initialize physics for the sake of MultiPhase routines
        sim().init_physics();

        // Initialize volume fraction field
        auto& vof = repo.get_field("vof");
        run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
            auto vof_arr = vof(lev).array(mfi);
            const auto& bx = mfi.validbox();
            initialize_volume_fractions(bx, vof_arr);
        });
        // Populate boundary cells
        vof.fillpatch(0.0);

        /* -- Check VOF boundary values from fillpatch -- */
        // Base level
        int lev = 0;
        int i = 0;
        int j = 0;
        int k = 0;
        for (amrex::MFIter mfi(vof(lev)); mfi.isValid(); ++mfi) {
            const auto& vof_arr = vof(lev).array(mfi);
            // Check lo and hi
            for (int i0 = 0; i0 < 2; ++i0) {
                int off = i0 * 3;
                // Small mesh, loop in serial for check
                for (int i1 = 0; i1 < 2; ++i1) {
                    for (int i2 = 0; i2 < 2; ++i2) {
                        // Organize indices
                        switch (dir) {
                        case 0:
                            i = -1 + off;
                            j = i1;
                            k = i2;
                            break;
                        case 1:
                            i = i1;
                            j = -1 + off;
                            k = i2;
                            break;
                        case 2:
                            i = i1;
                            j = i2;
                            k = -1 + off;
                            break;
                        }
                        // Calculate reference value
                        amrex::Real ref_val = m_vof_bdyval; // case 0
                        switch (vof_distr) {
                        case 1:
                            // hoextrap
                            ref_val = 1.0 - 0.1 * (i + j + k);
                            break;
                        case 2:
                            // foextrap
                            int ii = std::max(0, std::min(1, i));
                            int jj = std::max(0, std::min(1, j));
                            int kk = std::max(0, std::min(1, k));
                            ref_val = 1.0 - 0.1 * (ii + jj + kk);
                            break;
                        }
                        // Check against reference value
                        EXPECT_NEAR(vof_arr(i, j, k), ref_val, tol);
                    }
                }
            }
        }

        // Test positive and negative velocity
        for (int sign = -1; sign < 2; sign += 2) {

            // Get mac velocity fields and set values based on dir
            auto& umac = repo.get_field("u_mac");
            auto& vmac = repo.get_field("v_mac");
            auto& wmac = repo.get_field("w_mac");
            umac.setVal(dir == 0 ? m_vel * sign : 0.0);
            vmac.setVal(dir == 1 ? m_vel * sign : 0.0);
            wmac.setVal(dir == 2 ? m_vel * sign : 0.0);

            // Perform VOF solve
            // Advance states (new -> old)
            sim().pde_manager().advance_states();
            // Get equation handle and perform init
            auto& seqn = pde_mgr(
                amr_wind::pde::VOF::pde_name() + "-" +
                amr_wind::fvm::Godunov::scheme_name());
            seqn.initialize();
            seqn.compute_advection_term(amr_wind::FieldState::Old);
            seqn.post_solve_actions();

            // Get advected alpha
            const auto& advalpha_f = repo.get_field("advalpha_" + sdir);

            /* -- Check VOF boundary fluxes -- */
            for (amrex::MFIter mfi(vof(lev)); mfi.isValid(); ++mfi) {

                const auto& af = advalpha_f(lev).array(mfi);
                // Check lo and hi
                for (int i0 = 0; i0 < 2; ++i0) {
                    int off = i0 * 2;
                    // Small mesh, loop in serial for check
                    for (int i1 = 0; i1 < 2; ++i1) {
                        for (int i2 = 0; i2 < 2; ++i2) {
                            // Organize indices
                            switch (dir) {
                            case 0:
                                i = 0 + off;
                                j = i1;
                                k = i2;
                                break;
                            case 1:
                                i = i1;
                                j = 0 + off;
                                k = i2;
                                break;
                            case 2:
                                i = i1;
                                j = i2;
                                k = 0 + off;
                                break;
                            }
                            // Check whether flux is nonzero
                            amrex::Real advvof = 0.0;
                            amrex::Real advrho =
                                m_rho1 * advvof + m_rho2 * (1.0 - advvof);
                            if (nonzero_flux) {
                                EXPECT_GT(af(i, j, k), advrho);
                            } else {
                                EXPECT_EQ(af(i, j, k), advrho);
                            }
                        }
                    }
                }
            }
        }
    }
    const amrex::Real m_rho1 = 1000.0;
    const amrex::Real m_rho2 = 1.0;
    const amrex::Real m_vel = 5.0;
    const amrex::Real m_vof_bdyval = 1.0;
    const amrex::Real dt = 0.45 * 0.5 / m_vel; // first number is CFL
};

constexpr double tol1 = 1.0e-15;
constexpr double tol2 = 6.0e-2;

TEST_F(VOFBCTest, dirichletX) { testing_BC_coorddir(1, 0, tol1); }
TEST_F(VOFBCTest, slipwallY) { testing_BC_coorddir(2, 1, tol2); }
TEST_F(VOFBCTest, noslipwallZ) { testing_BC_coorddir(3, 2, tol1); }
TEST_F(VOFBCTest, pressureX) { testing_BC_coorddir(4, 0, tol1); }

} // namespace amr_wind_tests
