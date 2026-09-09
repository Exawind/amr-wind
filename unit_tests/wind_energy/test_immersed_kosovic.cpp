#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/test_utils.H"
#include "src/physics/ImmersedTerrain.H"
#include "src/turbulence/TurbulenceModel.H"
#include "src/utilities/math_ops.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

#include <limits>
#include <numbers>

using namespace amrex::literals;

namespace {
// 100 m plateau for x in [449, 576], flat ground elsewhere
void write_terrain(const std::string& fname)
{
    std::ofstream os(fname);
    os << "6\n2\n";
    os << "0.0\n448.0\n449.0\n576.0\n577.0\n1024.0\n";
    os << "0.0\n1024.0\n";
    os << "0.0\n0.0\n0.0\n0.0\n100.0\n100.0\n100.0\n100.0\n0.0\n0.0\n0.0\n0."
          "0\n";
}

void set_string(
    const std::string& prefix, const char* key, const std::string& v)
{
    amrex::ParmParse pp(prefix);
    pp.remove(key);
    pp.add(key, v);
}

//! Velocity (s z, v, 0) including ghost cells
void init_shear(
    kynema_sgf::Field& vel, const amrex::Real s, const amrex::Real v)
{
    const auto& mesh = vel.repo().mesh();
    const int nlevels = vel.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& varrs = vel(lev).arrays();
        amrex::ParallelFor(
            vel(lev), vel.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
                varrs[nbx](i, j, k, 0) = s * z;
                varrs[nbx](i, j, k, 1) = v;
                varrs[nbx](i, j, k, 2) = 0.0_rt;
            });
    }
    amrex::Gpu::streamSynchronize();
}
} // namespace

namespace kynema_sgf_tests {

/** Kosovic SGS viscosity with the ImmersedTerrain fields
 *
 *  Plateau terrain on a 32 x 32 x 16 mesh (dx = dy = dz = 32 m) with a
 *  shear flow (s z, v, 0). Checks the blanking inside the terrain, the
 *  log-law viscosity 2 rho u*^2 / |dU_t/dn| on the bottom and side wall
 *  patches, the (1 - beta) weight in a partial cell and the untouched
 *  Kosovic value away from the terrain.
 */
class ImmersedKosovicTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{32, 32, 16}};
            pp.addarr("n_cell", ncell);
            pp.add("blocking_factor", 2);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> probhi{{1024.0_rt, 1024.0_rt, 512.0_rt}};
            pp.addarr("prob_hi", probhi);
        }
        {
            amrex::ParmParse pp("turbulence");
            pp.add("model", (std::string) "Kosovic");
        }
        {
            amrex::ParmParse pp("incflo");
            amrex::Vector<std::string> physics{"ABL"};
            pp.addarr("physics", physics);
            pp.add("density", m_rho0);
            amrex::Vector<amrex::Real> vvec{8.0, 0.0, 0.0};
            pp.addarr("velocity", vvec);
            amrex::Vector<amrex::Real> gvec{0.0, 0.0, -9.81};
            pp.addarr("gravity", gvec);
        }
        {
            amrex::ParmParse pp("transport");
            pp.add("viscosity", 1.0e-5_rt);
            pp.add("reference_temperature", 300.0_rt);
        }
        {
            amrex::ParmParse pp("ABL");
            amrex::Vector<amrex::Real> t_hts{0.0, 100.0, 4000.0};
            pp.addarr("temperature_heights", t_hts);
            amrex::Vector<amrex::Real> t_vals{300.0, 300.0, 300.0};
            pp.addarr("temperature_values", t_vals);
            pp.add("surface_roughness_z0", m_z0);
        }
        {
            amrex::ParmParse pp("ImmersedTerrain");
            pp.add("uniform_roughness", m_z0);
        }
    }

    void setup()
    {
        write_terrain("terrain.amrwind");
        populate_parameters();
        initialize_mesh();
        auto& pde_mgr = sim().pde_manager();
        pde_mgr.register_icns();
        sim().init_physics();
        m_terrain = std::make_unique<Terrain>(sim());
        const int nlevels = sim().repo().num_active_levels();
        for (int lev = 0; lev < nlevels; ++lev) {
            m_terrain->initialize_fields(lev, sim().repo().mesh().Geom(lev));
        }
        sim().create_turbulence_model();
        sim().turbulence_model().post_init_actions();

        init_shear(sim().repo().get_field("velocity"), m_shear, m_vspan);
        sim().repo().get_field("density").setVal(m_rho0);
        sim().repo().get_field("temperature").setVal(300.0_rt);
        sim().turbulence_model().update_turbulent_viscosity(
            kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    }

    //! Kosovic value for the shear flow: rho Cs^2 Delta^2 |S| with |S| = s
    [[nodiscard]] amrex::Real kosovic_mu() const
    {
        const amrex::Real Cb = 0.36_rt;
        const amrex::Real Cs = std::sqrt(
            8.0_rt * (1.0_rt + Cb) /
            (27.0_rt * std::numbers::pi_v<amrex::Real> *
             std::numbers::pi_v<amrex::Real>));
        return m_rho0 * Cs * Cs * m_dz * m_dz * m_shear;
    }

    //! Neutral friction velocity from the reference speed at d2 = 1.5 dz
    [[nodiscard]] amrex::Real ustar(const amrex::Real m_ref) const
    {
        return 0.41_rt * m_ref / std::log(1.5_rt * m_dz / m_z0);
    }

    [[nodiscard]] amrex::Real speed(const int k) const
    {
        const amrex::Real z = (k + 0.5_rt) * m_dz;
        return std::sqrt((m_shear * z * m_shear * z) + (m_vspan * m_vspan));
    }

    using Terrain = kynema_sgf::immersedterrain::ImmersedTerrain;
    std::unique_ptr<Terrain> m_terrain;
    const amrex::Real m_dz{32.0_rt};
    const amrex::Real m_rho0{1.2_rt};
    const amrex::Real m_z0{0.1_rt};
    const amrex::Real m_shear{0.05_rt};
    const amrex::Real m_vspan{2.0_rt};
    const amrex::Real m_tol{
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt};
};

TEST_F(ImmersedKosovicTest, immersed_terrain_model)
{
    set_string("Kosovic", "terrain_model", "ImmersedTerrain");
    setup();
    const auto& mu = sim().repo().get_field("mu_turb");
    const auto& divNij = sim().repo().get_field("divNij");

    // Away from the terrain: plain Kosovic value
    EXPECT_NEAR(
        utils::field_probe(mu, 0, 5, 5, 8), kosovic_mu(), m_tol * kosovic_mu());
    // Inside the plateau: viscosity and non-linear term blanked
    EXPECT_NEAR(utils::field_probe(mu, 0, 15, 10, 1), 0.0_rt, m_tol);
    EXPECT_NEAR(utils::field_probe(divNij, 0, 15, 10, 1, 0), 0.0_rt, m_tol);

    // Partial cell on top of the plateau (beta = 0.125): one bottom patch,
    // reference cell two above, tangential speed gradient from the cells
    // above and below, weight 1 - beta
    {
        const amrex::Real us = ustar(speed(4));
        const amrex::Real dMdz = (speed(3) - speed(2)) / m_dz;
        const amrex::Real expected =
            (1.0_rt - 0.125_rt) * 2.0_rt * m_rho0 * us * us / dMdz;
        EXPECT_NEAR(
            utils::field_probe(mu, 0, 15, 10, 3), expected, m_tol * expected);
    }
    // Fluid cell beside the plateau wall: one x patch, the tangential speed
    // is the uniform v so the gradient sits at its floor of 0.01
    {
        const amrex::Real us = ustar(m_vspan);
        const amrex::Real expected = 2.0_rt * m_rho0 * us * us / 0.01_rt;
        EXPECT_NEAR(
            utils::field_probe(mu, 0, 13, 10, 1), expected, m_tol * expected);
    }
}

TEST_F(ImmersedKosovicTest, default_ignores_immersed_fields)
{
    // Default terrain_model = TerrainDrag: without terrain_blank the model
    // sees no terrain and the immersed fields are left alone
    setup();
    const auto& mu = sim().repo().get_field("mu_turb");
    EXPECT_NEAR(
        utils::field_probe(mu, 0, 5, 5, 8), kosovic_mu(), m_tol * kosovic_mu());
    EXPECT_NEAR(
        utils::field_probe(mu, 0, 15, 10, 1), kosovic_mu(),
        m_tol * kosovic_mu());
}

} // namespace kynema_sgf_tests
