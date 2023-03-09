#include "aw_test_utils/AmrexTest.H"
#include "amr-wind/incflo.H"

namespace amr_wind_tests {

namespace {

void init_density(amr_wind::Field& density, const int k_thresh = -3)
{
    const int nlevels = density.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {

        for (amrex::MFIter mfi(density(lev)); mfi.isValid(); ++mfi) {
            auto gbx = mfi.growntilebox();
            const auto& darr = density(lev).array(mfi);

            amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                darr(i, j, k) = (k > k_thresh) ? k + 4 : 0.0;
            });
        }
    }
}

amrex::Real get_Fgz_sum(amr_wind::Field& src_term)
{
    amrex::Real Fgz_sum = 0.0;

    for (int lev = 0; lev < src_term.repo().num_active_levels(); ++lev) {
        Fgz_sum += amrex::ReduceSum(
            src_term(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& Fg_arr) -> amrex::Real {
                amrex::Real Fgz_sum_fab = 0.0;

                amrex::Loop(
                    bx, [=, &Fgz_sum_fab](int i, int j, int k) noexcept {
                        Fgz_sum_fab += Fg_arr(i, j, k, 2);
                    });

                return Fgz_sum_fab;
            });
    }
    amrex::ParallelDescriptor::ReduceRealSum(Fgz_sum);
    return Fgz_sum;
}

void Fgtest_kernel(
    const amrex::Real Fz_ref,
    const int ncells,
    const int nz,
    const amr_wind::FieldState fstate,
    const bool make_ref_dens = false)
{
    incflo my_incflo;
    my_incflo.init_mesh();
    if (make_ref_dens) {
        // Create reference density field (only used if turned on via parser)
        auto& ref_dens =
            my_incflo.sim().repo().declare_field("reference_density", 1, 3, 1);
        init_density(ref_dens, nz / 2);
    }
    my_incflo.init_amr_wind_modules();
    auto& density = my_incflo.sim().repo().get_field("density").state(
        amr_wind::FieldState::NPH);
    auto& velocity = my_incflo.sim().repo().get_field("velocity");
    auto& grad_p = my_incflo.sim().repo().get_field("gp");
    auto& Fg_field = my_incflo.icns().fields().src_term;
    // Set density field
    init_density(density);
    // Set old density to unity, avoid NaN
    density.state(amr_wind::FieldState::Old).setVal(1.0);
    // Set zero velocity
    velocity.setVal(0.0);
    // Set zero pressure gradient
    grad_p.setVal(0.0);

    // Calculate forcing terms
    my_incflo.icns().compute_source_term(fstate);

    // Check forcing average value
    const amrex::Real Fz_avg = get_Fgz_sum(Fg_field) / ncells;
    EXPECT_NEAR(Fz_ref, Fz_avg, 1e-8);
}

} // namespace

class GravityForcingTest : public AmrexTest
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

            amrex::Vector<int> periodic{{1, 1, 1}};
            pp.addarr("is_periodic", periodic);
        }
        {
            amrex::ParmParse pp("incflo");
            pp.add("use_godunov", (int)1);
        }
        {
            amrex::ParmParse pp("ICNS");
            amrex::Vector<std::string> srcstr{"GravityForcing"};
            pp.addarr("source_terms", srcstr);
        }
    }

    const amrex::Real m_rho_0 = 2.0;
    const amrex::Real m_Fg = -9.81;
    const int m_nx = 2;
    const int m_ny = 2;
    const int m_nz = 16;
};

TEST_F(GravityForcingTest, full_term_u)
{
    // High-level setup
    populate_parameters();
    {
        amrex::ParmParse pp("ICNS");
        pp.add("use_perturb_pressure", (bool)false);
    }
    // Modify gravity to make sure it works
    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<amrex::Real> grav{0.0, 0.0, -5.0};
        pp.addarr("gravity", grav);
    }
    // Expected average gravity term
    amrex::Real Fg = -5.0;
    // Test with ordinary gravity term (rho not included)
    Fgtest_kernel(Fg, m_nx * m_ny * m_nz, m_nz, amr_wind::FieldState::Old);
}

TEST_F(GravityForcingTest, full_term_rhou)
{
    // High-level setup
    populate_parameters();
    {
        amrex::ParmParse pp("ICNS");
        pp.add("use_perturb_pressure", (bool)false);
    }
    // Modify gravity to make sure it works
    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<amrex::Real> grav{0.0, 0.0, -5.0};
        pp.addarr("gravity", grav);
    }
    // Expected average gravity term
    amrex::Real Fg = -5.0;
    amrex::Real fac = 0.0;
    for (int k = 0; k < m_nz; ++k) {
        fac += k + 4;
    }
    Fg *= fac / m_nz;
    // Test with ordinary gravity term (rho included)
    Fgtest_kernel(Fg, m_nx * m_ny * m_nz, m_nz, amr_wind::FieldState::New);
}

TEST_F(GravityForcingTest, perturb_const)
{
    const amrex::Real rho_ref = 0.5;
    const amrex::Real gz = -9.81;
    // High-level setup
    populate_parameters();
    {
        amrex::ParmParse pp("ICNS");
        pp.add("use_perturb_pressure", (bool)true);
    }
    // Modify gravity to make sure it works
    {
        amrex::ParmParse pp("incflo");
        pp.add("density", rho_ref);
    }
    // Expected average gravity term
    amrex::Real Fg = (1.0 - rho_ref) * gz / 1.0;
    // Test with ordinary gravity term (rho not multiplied)
    Fgtest_kernel(Fg, m_nx * m_ny * m_nz, m_nz, amr_wind::FieldState::Old);
}

TEST_F(GravityForcingTest, perturb_field)
{
    const amrex::Real gz = -9.81;
    // High-level setup
    populate_parameters();
    {
        amrex::ParmParse pp("ICNS");
        pp.add("use_perturb_pressure", (bool)true);
    }
    // Expected average gravity term
    amrex::Real Fg = gz;
    amrex::Real fac = 0.0;
    for (int k = 0; k < m_nz; ++k) {
        fac += (k > m_nz / 2) ? 0.0 : k + 4;
    }
    Fg *= fac / m_nz;
    // Test with ordinary gravity term (rho included)
    Fgtest_kernel(
        Fg, m_nx * m_ny * m_nz, m_nz, amr_wind::FieldState::New, true);
}

} // namespace amr_wind_tests