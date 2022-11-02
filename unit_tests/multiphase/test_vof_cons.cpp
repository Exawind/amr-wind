#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/equation_systems/vof/vof.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/utilities/tagging/CartBoxRefinement.H"

namespace amr_wind_tests {

//! Custom mesh class to directly specify refinement criteria
class NestRefineMesh : public AmrTestMesh
{
public:
    amrex::Vector<std::unique_ptr<amr_wind::RefinementCriteria>>&
    refine_criteria_vec()
    {
        return m_refine_crit;
    }

private:
    amrex::Vector<std::unique_ptr<amr_wind::RefinementCriteria>> m_refine_crit;
};

static void
initialize_volume_fractions(const int dir, const int nx, amr_wind::Field& vof)
{

    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        auto vof_arr = vof(lev).array(mfi);
        const auto& bx = mfi.validbox();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (dir < 0) {
                // Bottom left half is liquid, top right half is gas
                if (i + j + k == 3) {
                    vof_arr(i, j, k) = 0.5;
                } else {
                    if (i + j + k < 3) {
                        vof_arr(i, j, k) = 1.0;
                    } else {
                        vof_arr(i, j, k) = 0.0;
                    }
                }
            } else {
                int icheck = 0;
                // Left half is liquid, right half is gas
                switch (dir) {
                case 0:
                    icheck = i;
                    break;
                case 1:
                    icheck = j;
                    break;
                case 2:
                    icheck = k;
                    break;
                }
                if (2 * icheck + 1 == nx) {
                    vof_arr(i, j, k) = 0.5;
                } else {
                    if (2 * icheck + 1 < nx) {
                        vof_arr(i, j, k) = 1.0;
                    } else {
                        vof_arr(i, j, k) = 0.0;
                    }
                }
            }
        });
    });
    // Populate boundary cells
    vof.fillpatch(0.0);
}

static void initialize_adv_velocities(
    amr_wind::Field& vof,
    amr_wind::Field& umac,
    amr_wind::Field& vmac,
    amr_wind::Field& wmac,
    amrex::GpuArray<amrex::Real, 3> varr)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        auto um = umac(lev).array(mfi);
        auto vm = vmac(lev).array(mfi);
        auto wm = wmac(lev).array(mfi);
        const auto& gbx = mfi.growntilebox(1);
        amrex::ParallelFor(
            gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                um(i, j, k) = varr[0];
                vm(i, j, k) = varr[1];
                wm(i, j, k) = varr[2];
            });
    });
}

static void
check_accuracy(int dir, int nx, amrex::Real tol, amr_wind::Field& vof)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& vof_arr = vof(lev).const_array(mfi);

        // Loop manually through cells to check values
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < nx; ++j) {
                for (int k = 0; k < nx; ++k) {

                    int icheck = 0;
                    switch (dir) {
                    case 0:
                        icheck = i;
                        break;
                    case 1:
                        icheck = j;
                        break;
                    case 2:
                        icheck = k;
                        break;
                    }
                    // Check if current solution matches initial
                    // solution
                    if (2 * icheck + 1 == nx) {
                        EXPECT_NEAR(vof_arr(i, j, k), 0.5, tol);
                    } else {
                        if (2 * icheck + 1 < nx) {
                            EXPECT_NEAR(vof_arr(i, j, k), 1.0, tol);
                        } else {
                            EXPECT_NEAR(vof_arr(i, j, k), 0.0, tol);
                        }
                    }
                }
            }
        }
    });
}

class VOFConsTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{m_nx, m_nx, m_nx}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", m_nx);
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
        }
        {
            amrex::ParmParse pp("time");
            pp.add("fixed_dt", dt);
        }
        {
            amrex::ParmParse pp("VOF");
            pp.add("remove_debris", 0);
        }
    }

    void testing_coorddir(const int dir, amrex::Real CFL)
    {
        constexpr double tol = 1.0e-15;

        // Flow-through time
        const amrex::Real ft_time = 1.0 / m_vel;

        // Set timestep according to input
        dt = ft_time / ((amrex::Real)m_nx) * CFL;
        // Round to nearest integer timesteps
        int niter = (int)round(ft_time / dt);
        // Modify dt to fit niter
        dt = ft_time / ((amrex::Real)niter);

        populate_parameters();
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<int> periodic{{1, 1, 1}};
            pp.addarr("is_periodic", periodic);
        }
        // dir = -2 corresponds to multi-level
        if (dir == -2) {
            {
                amrex::ParmParse pp("amr");
                amrex::Vector<int> ncell{{m_nx + 1, m_nx + 1, m_nx + 1}};
                pp.add("max_level", 1);
                pp.add("max_grid_size", m_nx + 1);
                pp.add("blocking_factor", 2);
                pp.addarr("n_cell", ncell);
            }
            // Create the "input file"
            std::stringstream ss;
            ss << "1 // Number of levels" << std::endl;
            ss << "1 // Number of boxes at this level" << std::endl;
            ss << "0.8 0.5 0.5 0.9 0.5 0.5" << std::endl;

            create_mesh_instance<NestRefineMesh>();
            std::unique_ptr<amr_wind::CartBoxRefinement> box_refine(
                new amr_wind::CartBoxRefinement(sim()));
            box_refine->read_inputs(mesh(), ss);

            mesh<NestRefineMesh>()->refine_criteria_vec().push_back(
                std::move(box_refine));
        }

        initialize_mesh();

        auto& repo = sim().repo();

        // PDE manager, for access to VOF PDE later
        auto& pde_mgr = sim().pde_manager();
        // Setup of icns provides MAC velocities
        pde_mgr.register_icns();

        // Initialize physics for the sake of MultiPhase routines
        sim().init_physics();

        // Initialize volume fraction field
        auto& vof = repo.get_field("vof");
        initialize_volume_fractions(dir, m_nx, vof);
        // Get multiphase object
        auto& mphase = sim().physics_manager().get<amr_wind::MultiPhase>();

        // Initialize constant velocity field in single direction
        amrex::GpuArray<amrex::Real, 3> varr = {0};
        if (dir < 0) {
            varr[0] = m_vel;
            varr[1] = m_vel;
            varr[2] = m_vel;
        } else {
            varr[dir] = m_vel;
        }

        // Set advection velocities accordingly
        auto& umac = repo.get_field("u_mac");
        auto& vmac = repo.get_field("v_mac");
        auto& wmac = repo.get_field("w_mac");
        initialize_adv_velocities(vof, umac, vmac, wmac, varr);

        // Get initial VOF sum
        amrex::Real sum_vof0 = mphase.volume_fraction_sum();
        // Get equation handle and perform init
        auto& seqn = pde_mgr(
            amr_wind::pde::VOF::pde_name() + "-" +
            amr_wind::fvm::Godunov::scheme_name());
        seqn.initialize();

        for (int n = 0; n < niter; ++n) {
            // Perform VOF solve
            seqn.compute_advection_term(amr_wind::FieldState::Old);
            seqn.post_solve_actions();
            // Check conservation
            EXPECT_NEAR(mphase.volume_fraction_sum(), sum_vof0, tol);
        }

        if (dir >= 0) {
            check_accuracy(dir, m_nx, tol, vof);
        }
    }
    const amrex::Real m_rho1 = 1000.0;
    const amrex::Real m_rho2 = 1.0;
    const amrex::Real m_vel = 5.0;
    const int m_nx = 3;
    amrex::Real dt = 0.0; // will be set according to CFL
};

TEST_F(VOFConsTest, X) { testing_coorddir(0, 0.45); }
TEST_F(VOFConsTest, Y) { testing_coorddir(1, 0.45); }
TEST_F(VOFConsTest, Z) { testing_coorddir(2, 0.45); }
// Need multi-directional velocity and vof field to test communication of vof
// during directionally-split advection
TEST_F(VOFConsTest, CFL045) { testing_coorddir(-1, 0.45); }
TEST_F(VOFConsTest, CFL01) { testing_coorddir(-1, 0.1); }
// Test transport across multiple mesh levels - just check conservation
TEST_F(VOFConsTest, 2level) { testing_coorddir(-2, 0.5 * 0.45); }

} // namespace amr_wind_tests
