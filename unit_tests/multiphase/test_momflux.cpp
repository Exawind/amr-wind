#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"

namespace amr_wind_tests {

void init_field3(
    amr_wind::Field& fld,
    const amrex::Real in0,
    const amrex::Real in1,
    const amrex::Real in2)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();

    amrex::Real offset = 0.0;
    if (fld.field_location() == amr_wind::FieldLoc::CELL) offset = 0.5;

    for (int lev = 0; lev < nlevels; ++lev) {
        //const auto& dx = mesh.Geom(lev).CellSizeArray();
        //const auto& problo = mesh.Geom(lev).ProbLoArray();

        for (amrex::MFIter mfi(fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.growntilebox();
            const auto& farr = fld(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                //const amrex::Real x = problo[0] + (i + offset) * dx[0];
                //const amrex::Real y = problo[1] + (j + offset) * dx[1];
                //const amrex::Real z = problo[2] + (k + offset) * dx[2];

                farr(i, j, k, 0) = in0;
                farr(i, j, k, 1) = in1;
                farr(i, j, k, 2) = in2;
            });
        }
    }
}

class MassMomFluxOpTest : public MeshTest
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
            pp.add("use_godunov", (int) 1);
            pp.add("godunov_type", (std::string) "weno");
            //transport.model = TwoPhaseTransport
        }
        {
            amrex::ParmParse pp("time");
            pp.add("fixed_dt", dt);
        }
    }
    const amrex::Real m_rho1 = 1000.0;
    const amrex::Real m_rho2 = 1.0;
    const amrex::Real m_uvel = 5.0;
    const amrex::Real m_vvel = 0.0;
    const amrex::Real m_wvel = 0.0;
    const amrex::Real tol = 1e-8;
    const amrex::Real dt = 0.1 * 0.5 / m_uvel;
};

namespace {

void initialize_volume_fractions(
    const amrex::Geometry&,
    const amrex::Box& bx,
    amrex::Array4<amrex::Real>& vof_arr)
{

    // grow the box by 1 so that x,y,z go out of bounds and min(max()) corrects
    // it and it fills the ghosts with wall values
    amrex::ParallelFor(grow(bx, 1), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        if (i > 0) {
            vof_arr(i, j, k) = 0.0;
        } else {
            vof_arr(i, j, k) = 1.0;
        }
    });
    // Left half is liquid, right half is gas
}

} // namespace

TEST_F(MassMomFluxOpTest, fluxface)
{

    constexpr double tol = 1.0e-11;

    populate_parameters();
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<int> periodic{{1, 1, 1}};
        pp.addarr("is_periodic", periodic);
    }

    initialize_mesh();

    auto& repo = sim().repo();

    // Set up icns PDE, access to VOF PDE
    auto& pde_mgr = sim().pde_manager();
    auto& mom_eqn = pde_mgr.register_icns();
    mom_eqn.initialize();
    //vof_eqn.initialize();

    // Initialize physics for the sake of MultiPhase routines
    sim().init_physics();

    // Initialize constant velocity field
    auto& velocity = mom_eqn.fields().field.state(amr_wind::FieldState::Old);
    init_field3(velocity, m_uvel, m_vvel, m_wvel);

    // Initialize volume fraction field
    auto& vof = repo.get_field("vof");
    auto& geom = repo.mesh().Geom();
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        auto vof_arr = vof(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_volume_fractions(geom[lev], bx, vof_arr);
    });
    // Sync density field with vof
    auto& mphase = sim().physics_manager().get<amr_wind::MultiPhase>();
    mphase.set_density_via_vof();

    // Perform pre-advection step to get MAC velocity field (should be uniform)
    mom_eqn.pre_advection_actions(amr_wind::FieldState::Old);

    // Perform VOF solve
    for (auto& seqn : pde_mgr.scalar_eqns()) {
        if (seqn->fields().field.base_name() == "vof") {
            seqn->initialize();
            seqn->compute_advection_term(amr_wind::FieldState::Old);
            //seqn->post_solve_actions();
        }
    }

    // Perform momentum solve
    mom_eqn.compute_advection_term(amr_wind::FieldState::Old);

    // Get MAC velocities for use in testing
    auto& umac = repo.get_field("u_mac");
    auto& vmac = repo.get_field("v_mac");
    auto& wmac = repo.get_field("w_mac");

    // Get advected alpha
    auto& advalphax = repo.get_field("advalpha_x");

    // Get convective term
    auto& conv_term =
        mom_eqn.fields().conv_term.state(amr_wind::FieldState::New);

    // Base level
    int lev = 0;
    const auto& dx = geom[lev].CellSizeArray();
    const auto& problo = geom[lev].ProbLoArray();
    const auto& probhi = geom[lev].ProbHiArray();
    for (amrex::MFIter mfi(vof(lev)); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        amrex::Array4<amrex::Real const> const& u = umac(lev).array(mfi);
        amrex::Array4<amrex::Real const> const& v = vmac(lev).array(mfi);
        amrex::Array4<amrex::Real const> const& w = wmac(lev).array(mfi);
        amrex::Array4<amrex::Real const> const& afx = advalphax(lev).array(mfi);
        amrex::Array4<amrex::Real const> const& dqdt =
            conv_term(lev).array(mfi);

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // x is face location
                const amrex::Real x = problo[0] + i * dx[0];

                // UMAC
                EXPECT_NEAR(u(i, j, k), m_uvel, tol);

                // VMAC
                EXPECT_NEAR(v(i, j, k), m_vvel, tol);

                // WMAC
                EXPECT_NEAR(w(i, j, k), m_wvel, tol);

                // Test volume fractions at faces
                if (x == 0.5) {
                    // Center face (coming from left cell)
                    EXPECT_NEAR(afx(i, j, k), 1.0, tol);
                } else {
                    if (x == 0.0) {
                        // Left face (coming from right cell, periodic BC)
                        EXPECT_NEAR(afx(i, j, k), 0.0, tol);
                    }
                }

                // Test momentum fluxes by checking convective term
                if (i == 0) {
                    // Left cell (gas entering, liquid leaving)
                    EXPECT_NEAR(
                        dqdt(i, j, k, 0),
                        m_uvel * m_uvel * (m_rho2 - m_rho1) / 0.5, tol);
                } else {
                    if (i == 1) {
                        // Right cell (liquid entering, gas leaving)
                        EXPECT_NEAR(
                            dqdt(i, j, k, 0),
                            m_uvel * m_uvel * (m_rho1 - m_rho2) / 0.5, tol);
                    }
                }
            });
    }
}

} // namespace amr_wind_tests
