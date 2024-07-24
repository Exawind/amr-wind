
#include <utility>

#include "aw_test_utils/MeshTest.H"

#include "amr-wind/utilities/sampling/WaveEnergy.H"

namespace amr_wind_tests {

namespace {

void init_velocity(amr_wind::Field& fld)
{
    const int nlevels = fld.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& farrs = fld(lev).arrays();
        amrex::ParallelFor(
            fld(lev), amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                farrs[nbx](i, j, k, 0) = j;
                farrs[nbx](i, j, k, 1) = k;
                farrs[nbx](i, j, k, 2) = i;
            });
    }
    amrex::Gpu::synchronize();
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
        : amr_wind::wave_energy::WaveEnergy(sim, std::move(label))
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
            amrex::Vector<int> ncell{{m_nx, m_nx, m_nx}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", m_nx);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            pp.addarr("prob_lo", m_problo);
            pp.addarr("prob_hi", m_probhi);
        }
        {
            amrex::ParmParse pp("waveenergy");
            pp.add("output_frequency", 1);
            pp.add("water_level", m_wlev);
        }
        {
            amrex::ParmParse pp("incflo");
            amrex::Vector<amrex::Real> gvec{{0.0, 0.0, m_g}};
            pp.addarr("gravity", gvec);
            amrex::Vector<std::string> physics{"MultiPhase"};
            pp.addarr("physics", physics);
        }
        {
            amrex::ParmParse pp("MultiPhase");
            pp.add("density_fluid1", m_rho1);
        }
    }
    // Parameters
    const amrex::Vector<amrex::Real> m_problo{{0.0, 0.0, 0.0}};
    const amrex::Vector<amrex::Real> m_probhi{{2.0, 2.0, 1.0}};
    const int m_nx = 5;
    const amrex::Real m_wlev = 0.25;
    const amrex::Real m_g = -9;
    const amrex::Real m_rho1 = 888;
    const amrex::Real m_tol = 1e-12;
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
    const amrex::Real dx = 2.0 / (amrex::Real)m_nx;
    const amrex::Real dz = 1.0 / (amrex::Real)m_nx;
    const amrex::Real cell_vol = dx * dx * dz;
    amrex::Real ke_ref =
        0.5 * cell_vol / (m_wlev * 2.0 * 2.0) *
        (4.0 * 5.0 * (1.0 + 4.0 + 9.0 + 16.0) + 25.0 +
         0.5 * (4.0 * 15.0 + 5.0 * (4.0 + 16.0) +
                3.0 * (1.0 + 4.0 + 9.0 + 16.0)) +
         (4.0 * 10.0 + 5.0 * (1.0 + 9.0) + 2.0 * (1.0 + 4.0 + 9.0 + 16.0)));
    EXPECT_NEAR(ke, ke_ref, m_tol);
    // Formula has been integrated in z, and uses exact interface locations
    amrex::Real pe_exact =
        dx * dx * (-m_g) * 0.5 / (m_wlev * 2.0 * 2.0) *
            (15.0 * std::pow(2.5 * dz, 2) + 10.0 * std::pow(3.0 * dz, 2)) +
        0.5 * (-m_g) * m_wlev;
    EXPECT_NEAR(pe, pe_exact, m_tol);
}

} // namespace amr_wind_tests
