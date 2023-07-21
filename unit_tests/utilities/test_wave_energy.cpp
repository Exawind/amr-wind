
#include "aw_test_utils/MeshTest.H"

#include "amr-wind/utilities/sampling/WaveEnergy.H"

namespace amr_wind_tests {

namespace {

void init_velocity(amr_wind::Field& fld)
{
    const int nlevels = fld.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {

        for (amrex::MFIter mfi(fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.validbox();
            const auto& farr = fld(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                farr(i, j, k, 0) = j;
                farr(i, j, k, 1) = k;
                farr(i, j, k, 2) = i;
            });
        }
    }
}

void init_vof(amr_wind::Field& fld)
{
    const int nlevels = fld.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {

        for (amrex::MFIter mfi(fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.validbox();
            const auto& vof_arr = fld(lev).array(mfi);
            // Top half is air, bottom is gas, middle row alternates between
            // fully liquid and half liquid
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (k < 2) {
                    vof_arr(i, j, k) = 1.0;
                } else {
                    if (k == 2) {
                        if (i % 2 == 0) {
                            vof_arr(i, j, k) = 0.5;
                        } else {
                            vof_arr(i, j, k) = 1.0;
                        }
                    } else {
                        vof_arr(i, j, k) = 0.0;
                    }
                }
            });
        }
    }
}

class WaveEnergyImpl : public amr_wind::wave_energy::WaveEnergy
{
public:
    // cppcheck-suppress passedByValue
    WaveEnergyImpl(amr_wind::CFDSim& sim, std::string label)
        : amr_wind::wave_energy::WaveEnergy(sim, label)
    {}

protected:
    // No file output during test
    void prepare_ascii_file() override {}
    void write_ascii() override {}
};

} // namespace

class WaveEnergyTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{nx, nx, nx}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", nx);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
        {
            amrex::ParmParse pp("waveenergy");
            pp.add("output_frequency", 1);
            pp.add("water_level", wlev);
        }
        {
            amrex::ParmParse pp("incflo");
            amrex::Vector<amrex::Real> gvec{{0.0, 0.0, g}};
            pp.addarr("gravity", gvec);
            amrex::Vector<std::string> physics{"MultiPhase"};
            pp.addarr("physics", physics);
        }
        {
            amrex::ParmParse pp("MultiPhase");
            pp.add("density_fluid1", rho1);
        }
    }
    // Parameters
    const amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
    const amrex::Vector<amrex::Real> probhi{{2.0, 2.0, 1.0}};
    const int nx = 5;
    const amrex::Real wlev = 0.25;
    const amrex::Real g = -9;
    const amrex::Real rho1 = 888;
    const amrex::Real tol = 1e-12;
};

TEST_F(WaveEnergyTest, checkoutput)
{
    initialize_mesh();
    auto& repo = sim().repo();

    // Initialize physics for MultiPhase, icns must be registered too
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();

    // Get fields and initialize them
    auto& vof = repo.get_field("vof");
    auto& vel = repo.get_field("velocity");
    init_vof(vof);
    init_velocity(vel);

    // Initialize postprocessing tool
    WaveEnergyImpl tool(sim(), "waveenergy");
    tool.initialize();
    tool.post_advance_work();

    // Get answers
    amrex::Real ke = 0.0;
    amrex::Real pe = 0.0;
    tool.wave_energy(ke, pe);

    // Check answers
    const amrex::Real dx = 2.0 / (amrex::Real)nx;
    const amrex::Real dz = 1.0 / (amrex::Real)nx;
    const amrex::Real cell_vol = dx * dx * dz;
    amrex::Real ke_ref =
        0.5 * cell_vol / (wlev * 2.0 * 2.0) *
        (4.0 * 5.0 * (1.0 + 4.0 + 9.0 + 16.0) + 25.0 +
         0.5 * (4.0 * 15.0 + 5.0 * (4.0 + 16.0) +
                3.0 * (1.0 + 4.0 + 9.0 + 16.0)) +
         (4.0 * 10.0 + 5.0 * (1.0 + 9.0) + 2.0 * (1.0 + 4.0 + 9.0 + 16.0)));
    EXPECT_NEAR(ke, ke_ref, tol);
    // Formula has been integrated in z, and uses exact interface locations
    amrex::Real pe_exact =
        dx * dx * (-g) * 0.5 / (wlev * 2.0 * 2.0) *
            (15.0 * std::pow(2.5 * dz, 2) + 10.0 * std::pow(3.0 * dz, 2)) +
        0.5 * (-g) * wlev;
    EXPECT_NEAR(pe, pe_exact, tol);
}

} // namespace amr_wind_tests
