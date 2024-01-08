#include "aw_test_utils/MeshTest.H"
#include "amr-wind/physics/multiphase/hydrostatic_ops.H"

namespace amr_wind_tests {
namespace {

amrex::Real density_test_impl(
    amr_wind::MultiLevelVector& rho0,
    const amrex::Vector<amrex::Geometry> geom,
    const amrex::Real rho1,
    const amrex::Real rho2,
    const amrex::Real wlev)
{
    amrex::Real error_total = 0;

    for (int lev = 0; lev < rho0.size(); ++lev) {

        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();

        const auto ncells = rho0.ncells(lev);
        const auto* rho0_ptr = rho0.device_data(lev).data();

        error_total += amrex::Reduce::Sum(
            ncells,
            [=] AMREX_GPU_DEVICE(int k) {
                const amrex::Real zbtm = problo[2] + k * dx[2];
                amrex::Real vof = (wlev - zbtm) / dx[2];
                vof = amrex::max(vof, 0.0);
                vof = amrex::min(vof, 1.0);
                amrex::Real dens = vof * rho1 + (1.0 - vof) * rho2;
                return std::abs(rho0_ptr[k] - dens);
            },
            0.0);
    }
    return error_total;
}

amrex::Real pressure_test_impl(
    amr_wind::MultiLevelVector& p0,
    const amrex::Vector<amrex::Geometry> geom,
    const amrex::Real rho1,
    const amrex::Real rho2,
    const amrex::Real wlev,
    const amrex::Real gz)
{
    amrex::Real error_total = 0;

    for (int lev = 0; lev < p0.size(); ++lev) {

        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();
        const auto& probhi = geom[lev].ProbHiArray();

        const amrex::Real ht_max = probhi[2] - problo[2];
        const amrex::Real ht_min = 0.0;

        const auto ncells = p0.ncells(lev);
        const auto* p0_ptr = p0.device_data(lev).data();

        error_total += amrex::Reduce::Sum(
            ncells,
            [=] AMREX_GPU_DEVICE(int k) {
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
                const amrex::Real irhogz = -gz * (rho1 * ht_l + rho2 * ht_g);
                return std::abs(p0_ptr[k] - irhogz);
            },
            0.0);
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

    amr_wind::MultiLevelVector rho0{amr_wind::FieldLoc::CELL};
    rho0.resize(2, mesh().Geom());
    amr_wind::hydrostatic::define_rho0(
        rho0, m_rho1, m_rho2, m_wlev, mesh().Geom());

    amrex::Real error_total =
        density_test_impl(rho0, sim().mesh().Geom(), m_rho1, m_rho2, m_wlev);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, 1e-8);
}

TEST_F(MultiPhaseHydroStatic, reference_pressure)
{
    populate_parameters();
    initialize_mesh();

    amr_wind::MultiLevelVector p0{amr_wind::FieldLoc::NODE};
    p0.resize(2, mesh().Geom());
    amr_wind::hydrostatic::define_p0(
        p0, m_rho1, m_rho2, m_wlev, m_gz, mesh().Geom());

    amrex::Real error_total =
        pressure_test_impl(p0, mesh().Geom(), m_rho1, m_rho2, m_wlev, m_gz);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0, 1e-8);
}

} // namespace amr_wind_tests
