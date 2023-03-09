#include "aw_test_utils/AmrexTest.H"
#include "amr-wind/incflo.H"

namespace amr_wind_tests {

namespace {

void init_vel_z(amr_wind::Field& vel, const amrex::Real w_const)
{
    const int nlevels = vel.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {

        for (amrex::MFIter mfi(vel(lev)); mfi.isValid(); ++mfi) {
            auto gbx = mfi.growntilebox();
            const auto& varr = vel(lev).array(mfi);

            amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                varr(i, j, k, 2) = w_const;
            });
        }
    }
}

void init_ref_p(
    amr_wind::Field& ref_p,
    const amrex::Vector<amrex::Geometry> geom,
    const amrex::Real F_g,
    const amrex::Real rho_0)
{
    const int nlevels = ref_p.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = geom[lev].CellSizeArray();
        const auto& probhi = geom[lev].ProbHiArray();
        for (amrex::MFIter mfi(ref_p(lev)); mfi.isValid(); ++mfi) {
            auto nbx = mfi.nodaltilebox();
            const auto& p0_arr = ref_p(lev).array(mfi);

            amrex::ParallelFor(nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                // Height of pressure node
                const amrex::Real hnode = k * dx[2];
                // Integrated density from top
                const amrex::Real irho = rho_0 * (probhi[2] - hnode);

                // Multiply with force to get hydrostatic pressure
                p0_arr(i, j, k) = -irho * F_g;
            });
        }
    }
}

amrex::Real get_pbottom(amr_wind::Field& pressure)
{
    amrex::Real pb_sum = 0.0;

    for (int lev = 0; lev < pressure.repo().num_active_levels(); ++lev) {
        pb_sum += amrex::ReduceSum(
            pressure(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& nbx,
                amrex::Array4<amrex::Real const> const& p_arr) -> amrex::Real {
                amrex::Real pb_sum_fab = 0.0;

                amrex::Loop(
                    nbx, [=, &pb_sum_fab](int i, int j, int k) noexcept {
                        pb_sum_fab += (k == 0) ? p_arr(i, j, k) : 0.0;
                    });

                return pb_sum_fab;
            });
    }
    amrex::ParallelDescriptor::ReduceRealSum(pb_sum);
    return pb_sum;
}

void ptest_kernel(
    const amrex::Real rho_0,
    const amrex::Real w_0,
    const amrex::Real p_0,
    const int nbottom,
    const amrex::Real Fg = 0.0)
{
    incflo my_incflo;
    my_incflo.init_mesh();
    auto& density = my_incflo.sim().repo().get_field("density");
    auto& velocity = my_incflo.sim().repo().get_field("velocity");
    // Set uniform density
    density.setVal(rho_0);
    // Set velocity as it would be with gravity forcing
    init_vel_z(velocity, w_0);

    // If requested, form reference_pressure field
    if (Fg != 0.0) {
        // Pressure has 3 ghost points
        auto& p_ref_field = my_incflo.sim().repo().declare_nd_field(
            "reference_pressure", 1, 3, 1);
        init_ref_p(p_ref_field, my_incflo.sim().mesh().Geom(), Fg, rho_0);
    }

    // Time is set to non-zero: not testing initialization
    // Delta t is set to non-zero for result to work
    const amrex::Real time = 1.0;
    const amrex::Real dt = 1.0;
    // Apply projection
    my_incflo.ApplyProjection((density).vec_const_ptrs(), time, dt, false);
    // Get result
    auto& p = my_incflo.sim().repo().get_field("p");
    // Check result
    const amrex::Real pbottom = get_pbottom(p) / nbottom;
    EXPECT_NEAR(p_0, pbottom, 1e-8);
}

} // namespace

class ProjPerturb : public AmrexTest
{
protected:
    void populate_parameters()
    {
        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{m_nx, m_ny, m_nz}};
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
            amrex::ParmParse pp("incflo");
            pp.add("use_godunov", (int)1);
        }

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
        ppzhi.add("type", (std::string) "pressure_inflow");
    }

    const amrex::Real m_rho_0 = 1.0;
    const amrex::Real m_Fg = -9.81;
    const int m_nx = 2;
    const int m_ny = 2;
    const int m_nz = 16;
};

TEST_F(ProjPerturb, dynamic_only)
{
    // High-level setup
    populate_parameters();
    // Test with gravity term omitted
    ptest_kernel(m_rho_0, 0.0, 0.0, (m_nx + 1) * (m_ny + 1));
}

TEST_F(ProjPerturb, full_pressure)
{
    // High-level setup
    populate_parameters();
    // Test with gravity term included
    ptest_kernel(m_rho_0, m_Fg, -m_Fg, (m_nx + 1) * (m_ny + 1));
}

TEST_F(ProjPerturb, full_p_perturb)
{
    // High-level setup
    populate_parameters();
    // Test with gravity term omitted, then added as reference pressure
    ptest_kernel(m_rho_0, 0.0, -m_Fg, (m_nx + 1) * (m_ny + 1), m_Fg);
}

} // namespace amr_wind_tests