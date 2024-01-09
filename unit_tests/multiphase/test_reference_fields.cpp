#include "aw_test_utils/MeshTest.H"
#include "amr-wind/physics/multiphase/hydrostatic_ops.H"

namespace amr_wind_tests {
namespace {

amrex::Real density_test_impl(
    amr_wind::Field& rho0,
    const amrex::Vector<amrex::Geometry> geom,
    const amrex::Real rho1,
    const amrex::Real rho2,
    const amrex::Real wlev)
{
    amrex::Real error_total = 0;

    for (int lev = 0; lev < rho0.repo().num_active_levels(); ++lev) {

        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();

        error_total += amrex::ReduceSum(
            rho0(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& rho0_arr)
                -> amrex::Real {
                amrex::Real error = 0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                    const amrex::Real zbtm = problo[2] + k * dx[2];
                    amrex::Real vof = (wlev - zbtm) / dx[2];
                    vof = amrex::max(vof, 0.0);
                    vof = amrex::min(vof, 1.0);
                    amrex::Real dens = vof * rho1 + (1.0 - vof) * rho2;
                    error += std::abs(rho0_arr(i, j, k) - dens);
                });

                return error;
            });
    }
    return error_total;
}

amrex::Real pressure_test_impl(
    amr_wind::Field& p0,
    const amrex::Vector<amrex::Geometry> geom,
    const amrex::Real rho1,
    const amrex::Real rho2,
    const amrex::Real wlev,
    const amrex::Real gz,
    const int ngrow)
{
    amrex::Real error_total = 0;

    for (int lev = 0; lev < p0.repo().num_active_levels(); ++lev) {

        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();
        const auto& probhi = geom[lev].ProbHiArray();

        const amrex::Real ht_max = probhi[2] - problo[2];
        const amrex::Real ht_min = 0.0;

        error_total += amrex::ReduceSum(
            p0(lev), ngrow,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& nbx,
                amrex::Array4<amrex::Real const> const& p0_arr) -> amrex::Real {
                amrex::Real error = 0;

                amrex::Loop(nbx, [=, &error](int i, int j, int k) noexcept {
                    const amrex::Real znode = problo[2] + k * dx[2];
                    amrex::Real ht_g = probhi[2] - wlev;
                    amrex::Real ht_l = wlev - problo[2];
                    // Limit by location
                    ht_g = amrex::min(ht_g, probhi[2] - znode);
                    ht_l = amrex::min(ht_l, wlev - znode);
                    // Limit by bounds
                    ht_g = amrex::min(amrex::max(ht_g, ht_min), ht_max);
                    ht_l = amrex::min(amrex::max(ht_l, ht_min), ht_max);
                    // Integrated (-rho*g*z)
                    const amrex::Real irhogz =
                        -gz * (rho1 * ht_l + rho2 * ht_g);
                    error += std::abs(p0_arr(i, j, k) - irhogz);
                });

                return error;
            });
    }
    return error_total;
}

} // namespace

class MultiPhaseHydroStatic : public MeshTest
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
    }

    const amrex::Real m_rho1 = 1000.0;
    const amrex::Real m_rho2 = 1.0;
    const amrex::Real m_wlev = 0.5;
    const amrex::Real m_gz = -9.81;
    const int m_nx = 3;
};

TEST_F(MultiPhaseHydroStatic, reference_density)
{
    populate_parameters();
    initialize_mesh();
    auto& repo = sim().repo();
    const int ncomp = 1;
    const int nghost = 0;
    auto& rho0 = repo.declare_field("reference_density", ncomp, nghost);

    amr_wind::hydrostatic::define_rho0(
        rho0, m_rho1, m_rho2, m_wlev, sim().mesh().Geom());

    amrex::Real error_total =
        density_test_impl(rho0, sim().mesh().Geom(), m_rho1, m_rho2, m_wlev);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, 1e-8);
}

TEST_F(MultiPhaseHydroStatic, reference_pressure)
{
    populate_parameters();
    initialize_mesh();
    auto& repo = sim().repo();
    const int ncomp = 1;
    const int nghost = 3;
    auto& p0 = repo.declare_nd_field("reference_pressure", ncomp, nghost);

    amr_wind::hydrostatic::define_p0(
        p0, m_rho1, m_rho2, m_wlev, m_gz, sim().mesh().Geom());

    amrex::Real error_total = pressure_test_impl(
        p0, sim().mesh().Geom(), m_rho1, m_rho2, m_wlev, m_gz, nghost);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, 1e-8);
}

} // namespace amr_wind_tests
