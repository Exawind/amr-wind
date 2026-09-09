#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/test_utils.H"
#include "src/physics/ImmersedTerrain.H"
#include "src/turbulence/TurbulenceModel.H"
#include "src/equation_systems/tke/source_terms/KransAxell.H"
#include "src/utilities/constants.H"
#include "src/utilities/math_ops.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

#include <limits>

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

// z u v T tke, read by KransAxell for the mesoscale sponge
void write_rans_profile(const std::string& fname)
{
    std::ofstream os(fname);
    os << "0 8 0 300 0.1\n1000 8 0 300 0.1\n";
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

/** KLAxell eddy viscosity and KransAxell TKE source with the ImmersedTerrain
 *  fields: plateau terrain on a 32 x 32 x 16 mesh (dx = dy = dz = 32 m),
 *  neutral shear flow (s z, v, 0), uniform TKE.
 */
class ImmersedKLAxellTest : public MeshTest
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
            pp.add("model", (std::string) "KLAxell");
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
            pp.add("surface_temp_rate", 0.0_rt);
            pp.add("initial_wind_profile", true);
            pp.add("rans_1dprofile_file", (std::string) "rans_1d.info");
            amrex::Vector<amrex::Real> hts{0.0_rt, 100.0_rt, 4000.0_rt};
            pp.addarr("temperature_heights", hts);
            pp.addarr("wind_heights", hts);
            amrex::Vector<amrex::Real> t_vals{300.0_rt, 300.0_rt, 300.0_rt};
            pp.addarr("temperature_values", t_vals);
            amrex::Vector<amrex::Real> u_vals{8.0_rt, 8.0_rt, 8.0_rt};
            pp.addarr("u_values", u_vals);
            amrex::Vector<amrex::Real> v_vals{0.0_rt, 0.0_rt, 0.0_rt};
            pp.addarr("v_values", v_vals);
            amrex::Vector<amrex::Real> tke_vals{0.1_rt, 0.1_rt, 0.1_rt};
            pp.addarr("tke_values", tke_vals);
            pp.add("surface_temp_flux", 0.0_rt);
            pp.add("surface_roughness_z0", m_z0);
            // Keeps the neutral length scale and disables the sponge
            pp.add("meso_sponge_start", 1.0e5_rt);
        }
        {
            amrex::ParmParse pp("ImmersedTerrain");
            pp.add("uniform_roughness", m_z0);
        }
    }

    void setup()
    {
        write_terrain("terrain.amrwind");
        write_rans_profile("rans_1d.info");
        populate_parameters();
        initialize_mesh();
        auto& pde_mgr = sim().pde_manager();
        pde_mgr.register_icns();
        sim().init_physics();
        sim().create_transport_model();
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
        sim().repo().get_field("tke").setVal(m_tke);
        sim().repo().get_field("turb_lscale").setVal(10.0_rt);
        sim().time().delta_t() = m_dt;
        sim().turbulence_model().update_turbulent_viscosity(
            kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    }

    //! Neutral surface-layer mixing length at height z above the terrain
    [[nodiscard]] static amrex::Real lscale(const amrex::Real z)
    {
        const amrex::Real lambda = 30.0_rt;
        const amrex::Real kappa = 0.41_rt;
        return (lambda * kappa * z) / (lambda + (kappa * z));
    }

    //! Neutral KLAxell viscosity: rho Cmu l sqrt(k)
    [[nodiscard]] amrex::Real mu_rans(const amrex::Real z) const
    {
        return m_rho0 * m_Cmu * lscale(z) * std::sqrt(m_tke);
    }

    //! Dissipation Cmu^3 k^1.5 / (l + eps)
    [[nodiscard]] amrex::Real dissip(const amrex::Real z) const
    {
        return kynema_sgf::utils::powi(m_Cmu, 3) * std::pow(m_tke, 1.5_rt) /
               (lscale(z) + kynema_sgf::constants::EPS);
    }

    //! Neutral friction velocity from the reference speed at d2 = 1.5 dz
    [[nodiscard]] amrex::Real ustar(const amrex::Real m_ref) const
    {
        return 0.41_rt * m_ref / std::log(1.5_rt * m_dz / m_z0);
    }

    //! Log-law TKE target, u*^2 / Cmu^2 for zero heat flux
    [[nodiscard]] amrex::Real tke_exact(const amrex::Real m_ref) const
    {
        const amrex::Real us = ustar(m_ref);
        return us * us / (m_Cmu * m_Cmu);
    }

    [[nodiscard]] amrex::Real speed(const int k) const
    {
        const amrex::Real z = (k + 0.5_rt) * m_dz;
        return std::sqrt((m_shear * z * m_shear * z) + (m_vspan * m_vspan));
    }

    [[nodiscard]] amrex::Real drag_rate(const amrex::Real w_solid) const
    {
        const amrex::Real C = w_solid * 10.0_rt / m_dz;
        return (1.0_rt - std::exp(-C * m_dt)) / m_dt;
    }

    using Terrain = kynema_sgf::immersedterrain::ImmersedTerrain;
    std::unique_ptr<Terrain> m_terrain;
    const amrex::Real m_dz{32.0_rt};
    const amrex::Real m_dt{0.5_rt};
    const amrex::Real m_rho0{1.2_rt};
    const amrex::Real m_z0{0.1_rt};
    const amrex::Real m_shear{0.05_rt};
    const amrex::Real m_vspan{2.0_rt};
    const amrex::Real m_tke{0.1_rt};
    const amrex::Real m_Cmu{0.556_rt};
    const amrex::Real m_tol{
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt};
};

TEST_F(ImmersedKLAxellTest, viscosity_immersed_terrain_model)
{
    set_string("KLAxell", "terrain_model", "ImmersedTerrain");
    setup();
    const auto& mu = sim().repo().get_field("mu_turb");

    // Away from the terrain: mixing length from the height above ground
    EXPECT_NEAR(
        utils::field_probe(mu, 0, 5, 5, 8), mu_rans(272.0_rt),
        m_tol * mu_rans(272.0_rt));
    // Inside the plateau: blanked
    EXPECT_NEAR(utils::field_probe(mu, 0, 15, 10, 1), 0.0_rt, m_tol);
    // Partial cell on top of the plateau (beta = 0.125): height above the
    // terrain floored at dz/2, weight 1 - beta
    {
        const amrex::Real expected = (1.0_rt - 0.125_rt) * mu_rans(16.0_rt);
        EXPECT_NEAR(
            utils::field_probe(mu, 0, 15, 10, 3), expected, m_tol * expected);
    }
}

TEST_F(ImmersedKLAxellTest, tke_source_immersed_terrain_model)
{
    set_string("KLAxell", "terrain_model", "ImmersedTerrain");
    setup();
    auto& src_term = sim().repo().get_field("tke_src_term");
    src_term.setVal(0.0_rt);
    kynema_sgf::pde::tke::KransAxell src(sim());
    src(0, kynema_sgf::FieldState::New, src_term(0));
    const amrex::Real tau = 5.0_rt * m_dt;

    // Inside the plateau: production and dissipation blanked, TKE damped at
    // the exact-integration rate of 10/dz
    EXPECT_NEAR(
        utils::field_probe(src_term, 0, 15, 10, 1), -drag_rate(1.0_rt) * m_tke,
        m_tol);
    // Partial cell on top of the plateau: (1 - beta) times production minus
    // dissipation, log-law relaxation from the cell two above, damping at
    // beta times the drag rate
    {
        const amrex::Real w = 0.125_rt;
        const amrex::Real mu_c = (1.0_rt - w) * mu_rans(16.0_rt);
        const amrex::Real main = (m_shear * m_shear * mu_c) - dissip(16.0_rt);
        const amrex::Real wall =
            (1.0_rt - w) * (tke_exact(speed(4)) - m_tke) / tau;
        const amrex::Real expected =
            ((1.0_rt - w) * main) + wall - (drag_rate(w) * m_tke);
        EXPECT_NEAR(
            utils::field_probe(src_term, 0, 15, 10, 3), expected,
            m_tol * std::abs(expected));
    }
    // Fluid cell on the bottom domain face: plain bottom wall function
    {
        const amrex::Real main =
            (m_shear * m_shear * mu_rans(16.0_rt)) - dissip(16.0_rt);
        const amrex::Real wall = (tke_exact(speed(1)) - m_tke) / tau;
        const amrex::Real expected = main + wall;
        EXPECT_NEAR(
            utils::field_probe(src_term, 0, 5, 5, 0), expected,
            m_tol * std::abs(expected));
    }
    // Fluid cell away from the walls: production minus dissipation only
    {
        const amrex::Real expected =
            (m_shear * m_shear * mu_rans(272.0_rt)) - dissip(272.0_rt);
        EXPECT_NEAR(
            utils::field_probe(src_term, 0, 5, 5, 8), expected,
            m_tol * std::abs(expected));
    }
}

TEST_F(ImmersedKLAxellTest, default_ignores_immersed_fields)
{
    // Default terrain_model = TerrainDrag: without terrain_blank the model
    // sees no terrain, so the solid cell keeps the flat-ground value
    setup();
    const auto& mu = sim().repo().get_field("mu_turb");
    EXPECT_NEAR(
        utils::field_probe(mu, 0, 15, 10, 1), mu_rans(48.0_rt),
        m_tol * mu_rans(48.0_rt));
}

} // namespace kynema_sgf_tests
