#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/iter_tools.H"
#include "ks_test_utils/test_utils.H"
#include "src/physics/ImmersedTerrain.H"
#include "src/physics/ImmersedWallModel.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

#include <array>
#include <cmath>
#include <limits>

using namespace amrex::literals;

/** Manufactured-solution test of the immersed wall model.
 *
 *  An exact neutral log-law velocity, tangential to the terrain and a
 *  function of the true (closest-point) distance to it, is imposed above a
 *  plane slope and above a Gaussian ridge. The wall-model geometry of
 *  ImmersedWallModel.H is then asked to recover the friction velocity and the
 *  log-law target in every surface cell. The relative errors are measured on
 *  three meshes for each wall-distance method. On the plane the
 *  surface_normal method with actual reference distances must be exact; on
 *  the ridge it must converge, while the face-based methods carry a
 *  slope-dependent error that does not vanish with the mesh.
 */
namespace {

constexpr amrex::Real kappa = 0.41_rt;
constexpr amrex::Real z0 = 0.1_rt;
constexpr amrex::Real ustar_exact = 0.5_rt;
constexpr amrex::Real Lx = 1024.0_rt;

//! Terrain shapes: plane slope and Gaussian ridge, uniform in y
struct Terrain
{
    int kind{0}; //!< 0 plane, 1 ridge
    amrex::Real h0{100.0_rt};
    amrex::Real slope{0.36397_rt}; //!< tan(20 deg)
    amrex::Real amp{200.0_rt};
    amrex::Real xc{512.0_rt};
    amrex::Real sig{200.0_rt};
    amrex::Real lx{Lx};

    [[nodiscard]] AMREX_GPU_HOST_DEVICE amrex::Real h(const amrex::Real x) const
    {
        if (kind == 0) {
            return h0 + (slope * x);
        }
        const amrex::Real a = (x - xc) / sig;
        return h0 + (amp * std::exp(-0.5_rt * a * a));
    }
    [[nodiscard]] AMREX_GPU_HOST_DEVICE amrex::Real
    dh(const amrex::Real x) const
    {
        if (kind == 0) {
            return slope;
        }
        const amrex::Real a = (x - xc) / sig;
        return -amp * a / sig * std::exp(-0.5_rt * a * a);
    }
    //! Closest point on the curve z = h(s) to (x, z): golden section search
    [[nodiscard]] AMREX_GPU_HOST_DEVICE amrex::Real
    closest_s(const amrex::Real x, const amrex::Real z) const
    {
        amrex::Real a = amrex::max<amrex::Real>(x - 600.0_rt, -200.0_rt);
        amrex::Real b = amrex::min<amrex::Real>(x + 600.0_rt, lx + 200.0_rt);
        auto dist2 = [&](const amrex::Real s) {
            const amrex::Real dz = z - h(s);
            return ((x - s) * (x - s)) + (dz * dz);
        };
        const amrex::Real gr = 0.6180339887498949_rt;
        amrex::Real c = b - (gr * (b - a));
        amrex::Real d = a + (gr * (b - a));
        for (int it = 0; it < 80; ++it) {
            if (dist2(c) < dist2(d)) {
                b = d;
            } else {
                a = c;
            }
            c = b - (gr * (b - a));
            d = a + (gr * (b - a));
        }
        return 0.5_rt * (a + b);
    }
};

void write_terrain(const std::string& fname, const Terrain& t)
{
    std::ofstream os(fname);
    os.precision(17);
    const int nx = 4097;
    os << nx << "\n2\n";
    for (int i = 0; i < nx; ++i) {
        os << Lx * i / (nx - 1) << "\n";
    }
    os << "0.0\n256.0\n";
    for (int i = 0; i < nx; ++i) {
        const amrex::Real hv = t.h(Lx * i / (nx - 1));
        os << hv << "\n" << hv << "\n";
    }
}

struct Config
{
    const char* label;
    kynema_sgf::immersed_wall::WallModel model;
    bool actual_reference;
    bool center_weight;
};

struct Errors
{
    amrex::Real ustar{0.0_rt};
    amrex::Real target{0.0_rt};
    amrex::Real count{0.0_rt};
};

} // namespace

namespace kynema_sgf_tests {

class ImmersedWallMMSTest : public MeshTest
{
public:
    void populate_parameters() override { MeshTest::populate_parameters(); }

    void setup(const int n, const Terrain& terrain)
    {
        write_terrain("terrain.amrwind", terrain);
        // Base defaults first so that the mesh parameters below take effect
        populate_parameters();
        {
            amrex::ParmParse pp("amr");
            pp.remove("n_cell");
            pp.remove("blocking_factor");
            pp.remove("max_grid_size");
            amrex::Vector<int> ncell{{n, 8, n / 2}};
            pp.addarr("n_cell", ncell);
            pp.add("blocking_factor", 2);
            pp.add("max_grid_size", 64);
        }
        {
            amrex::ParmParse pp("geometry");
            pp.remove("prob_hi");
            amrex::Vector<amrex::Real> probhi{{Lx, 256.0_rt, 512.0_rt}};
            pp.addarr("prob_hi", probhi);
        }
        initialize_mesh();
        sim().pde_manager().register_icns();
        sim().init_physics();
        m_terrain = std::make_unique<ImmersedTerrainT>(sim());
        m_terrain->initialize_fields(0, sim().repo().mesh().Geom(0));

        // Exact log-law velocity tangential to the terrain, including ghosts
        auto& vel = sim().repo().get_field("velocity");
        const auto& geom = sim().repo().mesh().Geom(0);
        const auto dx = geom.CellSizeArray();
        const auto plo = geom.ProbLoArray();
        const auto vel_arrs = vel(0).arrays();
        const Terrain tt = terrain;
        const amrex::Real z0_d = z0;
        const amrex::Real kappa_d = kappa;
        const amrex::Real ustar_d = ustar_exact;
        amrex::ParallelFor(
            vel(0), amrex::IntVect(1),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real x = plo[0] + ((i + 0.5_rt) * dx[0]);
                const amrex::Real z = plo[2] + ((k + 0.5_rt) * dx[2]);
                const amrex::Real s = tt.closest_s(x, z);
                const amrex::Real hs = tt.h(s);
                const amrex::Real d =
                    std::sqrt(((x - s) * (x - s)) + ((z - hs) * (z - hs)));
                const amrex::Real slope = tt.dh(s);
                const amrex::Real tn = std::sqrt(1.0_rt + (slope * slope));
                amrex::Real umag = 0.0_rt;
                if (z > hs && d > z0_d) {
                    umag = ustar_d / kappa_d * std::log(d / z0_d);
                }
                vel_arrs[nbx](i, j, k, 0) = umag / tn;
                vel_arrs[nbx](i, j, k, 1) = 0.0_rt;
                vel_arrs[nbx](i, j, k, 2) = umag * slope / tn;
            });
        amrex::Gpu::streamSynchronize();
    }

    Errors evaluate(const Config& cfg)
    {
        using namespace kynema_sgf::immersed_wall;
        const auto& repo = sim().repo();
        const auto& vel = repo.get_field("velocity")(0);
        const auto& frac = repo.get_field("terrain_fraction")(0);
        const auto& surf = repo.get_field("terrain_surface")(0);
        const auto& geom = repo.mesh().Geom(0);
        const auto dx = geom.CellSizeArray();
        const auto plo = geom.ProbLoArray();
        const amrex::Real thr = 0.5_rt;
        WallParams wp{};
        wp.kappa = kappa;
        const auto model = cfg.model;
        const amrex::Real z0_d = z0;
        const amrex::Real kappa_d = kappa;
        const amrex::Real ustar_d = ustar_exact;
        const bool actual = cfg.actual_reference;
        const bool center = cfg.center_weight;

        // 0: ustar error, 1: target error, 2: count
        Errors e;
        for (int which = 0; which < 3; ++which) {
            e_val(which) = amrex::ReduceSum(
                vel, frac, surf, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& v,
                    amrex::Array4<amrex::Real const> const& f,
                    amrex::Array4<amrex::Real const> const& s) -> amrex::Real {
                    amrex::Real acc = 0.0_rt;
                    amrex::Loop(bx, [=, &acc](int i, int j, int k) {
                        const amrex::Real beta = f(i, j, k, 0);
                        const bool surface =
                            (beta > 0.0_rt && beta < 1.0_rt) ||
                            (beta == 0.0_rt && f(i, j, k - 1, 0) >= thr);
                        if (!surface) {
                            return;
                        }
                        if (solid_weight(beta, center, thr) >= 1.0_rt) {
                            return;
                        }
                        const amrex::Real z_c = plo[2] + ((k + 0.5_rt) * dx[2]);
                        amrex::GpuArray<WallPatch, 2 * AMREX_SPACEDIM>
                            patches{};
                        const int np = wall_patches(
                            model, i, j, k, beta, f, s, dx, z_c, z0_d, thr,
                            actual, patches.data());
                        amrex::Real eu = 0.0_rt;
                        amrex::Real et = 0.0_rt;
                        amrex::Real wsum = 0.0_rt;
                        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> uc{
                            v(i, j, k, 0), v(i, j, k, 1), v(i, j, k, 2)};
                        const amrex::Real u_exact = magnitude(uc);
                        for (int ip = 0; ip < np; ++ip) {
                            const WallPatch& p = patches[ip];
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>
                                ur{v(p.ir, p.jr, p.kr, 0),
                                   v(p.ir, p.jr, p.kr, 1),
                                   v(p.ir, p.jr, p.kr, 2)};
                            const amrex::Real m_ref =
                                magnitude(tangential(ur, p.nrm));
                            const amrex::Real us =
                                friction_velocity(m_ref, p, z0_d, wp);
                            const amrex::Real dd1 =
                                amrex::max<amrex::Real>(p.d1, z0_d);
                            const amrex::Real target =
                                us / kappa_d * wp.phi_m(dd1, z0_d);
                            eu += p.weight * std::abs(us - ustar_d) / ustar_d;
                            et +=
                                p.weight * std::abs(target - u_exact) / ustar_d;
                            wsum += p.weight;
                        }
                        if (wsum > 0.0_rt) {
                            acc += (which == 0)   ? eu / wsum
                                   : (which == 1) ? et / wsum
                                                  : 1.0_rt;
                        }
                    });
                    return acc;
                });
            amrex::ParallelDescriptor::ReduceRealSum(e_val(which));
        }
        e.ustar = m_vals[0] / m_vals[2];
        e.target = m_vals[1] / m_vals[2];
        e.count = m_vals[2];
        return e;
    }

    amrex::Real& e_val(const int which) { return m_vals[which]; }

    using ImmersedTerrainT = kynema_sgf::immersedterrain::ImmersedTerrain;
    std::unique_ptr<ImmersedTerrainT> m_terrain;
    std::array<amrex::Real, 3> m_vals{{0.0_rt, 0.0_rt, 0.0_rt}};
};

namespace {
using kynema_sgf::immersed_wall::WallModel;
const std::array<Config, 5> configs{
    {{.label = "cell_offset",
      .model = WallModel::cell_offset,
      .actual_reference = false,
      .center_weight = false},
     {.label = "terrain_height",
      .model = WallModel::terrain_height,
      .actual_reference = false,
      .center_weight = false},
     {.label = "surface_normal (nominal d2)",
      .model = WallModel::surface_normal,
      .actual_reference = false,
      .center_weight = false},
     {.label = "surface_normal (actual d2)",
      .model = WallModel::surface_normal,
      .actual_reference = true,
      .center_weight = false},
     {.label = "surface_normal (actual d2, center)",
      .model = WallModel::surface_normal,
      .actual_reference = true,
      .center_weight = true}}};
constexpr std::array<int, 3> meshes{{32, 64, 128}};
} // namespace

namespace {
// errors[kind][config][mesh][ustar|target]; filled by the per-resolution
// tests below, which gtest runs in declaration order
std::array<std::array<std::array<std::array<amrex::Real, 2>, 3>, 5>, 2> g_err{};

void run_case(ImmersedWallMMSTest& t, const int kind, const int im)
{
    Terrain terrain;
    terrain.kind = kind;
    if (kind == 1) {
        terrain.h0 = 50.0_rt;
    }
    t.setup(meshes[im], terrain);
    for (int ic = 0; ic < 5; ++ic) {
        const Errors e = t.evaluate(configs[ic]);
        g_err[kind][ic][im][0] = e.ustar;
        g_err[kind][ic][im][1] = e.target;
        amrex::Print() << "MMS " << (kind == 0 ? "plane" : "ridge")
                       << "  n=" << meshes[im] << "  " << configs[ic].label
                       << "  ustar err=" << e.ustar
                       << "  target err=" << e.target << "  cells=" << e.count
                       << "\n";
    }
}
} // namespace

TEST_F(ImmersedWallMMSTest, plane_n32) { run_case(*this, 0, 0); }
TEST_F(ImmersedWallMMSTest, plane_n64) { run_case(*this, 0, 1); }
TEST_F(ImmersedWallMMSTest, plane_n128_orders)
{
    run_case(*this, 0, 2);
    const amrex::Real tol_exact =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e5_rt;
    // consistent normal method: exact on a plane (both drag weights)
    EXPECT_LT(g_err[0][3][2][0], tol_exact);
    EXPECT_LT(g_err[0][3][2][1], tol_exact);
    EXPECT_LT(g_err[0][4][2][0], tol_exact);
    // face-based and nominal methods: slope error does not converge
    for (int ic : {0, 1, 2}) {
        EXPECT_GT(g_err[0][ic][2][0], 0.5_rt * g_err[0][ic][0][0]);
        EXPECT_GT(g_err[0][ic][2][0], 1.0e-3_rt);
    }
}

TEST_F(ImmersedWallMMSTest, ridge_n32) { run_case(*this, 1, 0); }
TEST_F(ImmersedWallMMSTest, ridge_n64) { run_case(*this, 1, 1); }
TEST_F(ImmersedWallMMSTest, ridge_n128_orders)
{
    run_case(*this, 1, 2);
    for (int ic = 0; ic < 5; ++ic) {
        const amrex::Real order =
            std::log2(g_err[1][ic][1][0] / g_err[1][ic][2][0]);
        amrex::Print() << "MMS ridge  " << configs[ic].label
                       << "  ustar order (64->128) = " << order << "\n";
    }
    // The consistent normal method converges. Its per-cell error is the
    // tangent-plane approximation of the distance on a curved surface, which
    // scales like dx / ln(dx/z0), so the measured order is somewhat below 1
    // (0.77 over 32 -> 128 here). The nominal and face-based methods carry a
    // slope error independent of the mesh.
    const amrex::Real order_actual =
        0.5_rt * std::log2(g_err[1][3][0][0] / g_err[1][3][2][0]);
    EXPECT_GT(order_actual, 0.6_rt);
    EXPECT_LT(g_err[1][3][2][0], 0.02_rt * g_err[1][2][2][0]);
    EXPECT_LT(g_err[1][3][2][0], 0.002_rt * g_err[1][0][2][0]);
    for (int ic : {0, 1, 2}) {
        EXPECT_GT(g_err[1][ic][2][0], 0.5_rt * g_err[1][ic][0][0]);
    }
}

} // namespace kynema_sgf_tests
