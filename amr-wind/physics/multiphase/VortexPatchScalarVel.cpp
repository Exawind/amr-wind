#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/physics/multiphase/VortexPatchScalarVel.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

VortexPatchScalarVel::VortexPatchScalarVel(CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
    , m_density(sim.repo().get_field("density"))
{
    amrex::ParmParse pp(identifier());
    pp.queryarr("location", m_loc, 0, AMREX_SPACEDIM);
    pp.query("radius", m_radius);
    pp.query("smooth_factor", m_sfactor);
    pp.query("period", m_TT);
    amrex::ParmParse pinc("incflo");
    pinc.add("prescribe_velocity", true);
}

/** Initialize the velocity and levelset fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::VortexPatchScalarVelFieldInit
 */
void VortexPatchScalarVel::initialize_fields(
    int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& levelset = m_levelset(level);
    auto& density = m_density(level);
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();

    const amrex::Real xc = m_loc[0];
    const amrex::Real yc = m_loc[1];
    const amrex::Real zc = m_loc[2];
    const amrex::Real radius = m_radius;

    const auto& mphase = m_sim.physics_manager().get<MultiPhase>();
    const amrex::Real rho1 = mphase.rho1();
    const amrex::Real rho2 = mphase.rho2();

    auto& u_mac = m_sim.repo().get_field("u_mac")(level);
    auto& v_mac = m_sim.repo().get_field("v_mac")(level);
    auto& w_mac = m_sim.repo().get_field("w_mac")(level);

    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        auto uf = u_mac.array(mfi);
        auto vf = v_mac.array(mfi);
        auto wf = w_mac.array(mfi);
        auto vel = velocity.array(mfi);
        auto phi = levelset.array(mfi);
        auto rho = density.array(mfi);
        const amrex::Real eps = std::cbrt(2. * dx[0] * dx[1] * dx[2]);
        const amrex::Real eps_vel = radius * m_sfactor;

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                const amrex::Real xf = problo[0] + i * dx[0];
                const amrex::Real yf = problo[1] + j * dx[1];
                const amrex::Real zf = problo[2] + k * dx[2];
                uf(i, j, k) = 2.0 * std::sin(M_PI * xf) * std::sin(M_PI * xf) *
                              std::sin(2.0 * M_PI * y) *
                              std::sin(2.0 * M_PI * z);
                vf(i, j, k) = -std::sin(M_PI * yf) * std::sin(M_PI * yf) *
                              std::sin(2.0 * M_PI * x) *
                              std::sin(2.0 * M_PI * z);
                wf(i, j, k) = -std::sin(M_PI * zf) * std::sin(M_PI * zf) *
                              std::sin(2.0 * M_PI * x) *
                              std::sin(2.0 * M_PI * y);
                // Only the x component is nonzero
                vel(i, j, k, 1) = 0.0;
                vel(i, j, k, 2) = 0.0;

                phi(i, j, k) =
                    radius - std::sqrt(
                                 (x - xc) * (x - xc) + (y - yc) * (y - yc) +
                                 (z - zc) * (z - zc));
                amrex::Real smooth_heaviside;
                if (phi(i, j, k) > eps) {
                    smooth_heaviside = 1.0;
                } else if (phi(i, j, k) < -eps) {
                    smooth_heaviside = 0.;
                } else {
                    smooth_heaviside =
                        0.5 *
                        (1.0 + phi(i, j, k) / eps +
                         1.0 / M_PI * std::sin(phi(i, j, k) * M_PI / eps));
                }
                rho(i, j, k) =
                    rho1 * smooth_heaviside + rho2 * (1.0 - smooth_heaviside);

                // Smoothed step function for u velocity, which is treated as a
                // scalar, much more smoothed than levelset and not based on the
                // mesh size
                if (phi(i, j, k) > eps_vel) {
                    vel(i, j, k, 0) = 1.0;
                } else if (phi(i, j, k) < -eps_vel) {
                    vel(i, j, k, 0) = 0.;
                } else {
                    vel(i, j, k, 0) =
                        0.5 *
                        (1.0 + phi(i, j, k) / eps_vel +
                         1.0 / M_PI * std::sin(phi(i, j, k) * M_PI / eps_vel));
                }
            });
    }
    m_velocity.fillpatch(0.0);
}

void VortexPatchScalarVel::pre_advance_work()
{
    const auto& time =
        m_sim.time().current_time() + 0.5 * m_sim.time().deltaT();

    const int nlevels = m_sim.repo().num_active_levels();
    const auto& geom = m_sim.mesh().Geom();

    // Overriding the velocity field
    for (int lev = 0; lev < nlevels; ++lev) {
        auto& u_mac = m_sim.repo().get_field("u_mac")(lev);
        auto& v_mac = m_sim.repo().get_field("v_mac")(lev);
        auto& w_mac = m_sim.repo().get_field("w_mac")(lev);
        for (amrex::MFIter mfi(m_velocity(lev)); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.growntilebox(1);
            const auto& dx = geom[lev].CellSizeArray();
            const auto& problo = geom[lev].ProbLoArray();
            const amrex::Real TT = m_TT;
            auto uf = u_mac.array(mfi);
            auto vf = v_mac.array(mfi);
            auto wf = w_mac.array(mfi);
            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                    const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                    const amrex::Real xf = problo[0] + i * dx[0];
                    const amrex::Real yf = problo[1] + j * dx[1];
                    const amrex::Real zf = problo[2] + k * dx[2];
                    uf(i, j, k) =
                        2.0 * std::sin(M_PI * xf) * std::sin(M_PI * xf) *
                        std::sin(2.0 * M_PI * y) * std::sin(2.0 * M_PI * z) *
                        std::cos(M_PI * time / TT);
                    vf(i, j, k) = -std::sin(M_PI * yf) * std::sin(M_PI * yf) *
                                  std::sin(2.0 * M_PI * x) *
                                  std::sin(2.0 * M_PI * z) *
                                  std::cos(M_PI * time / TT);
                    wf(i, j, k) = -std::sin(M_PI * zf) * std::sin(M_PI * zf) *
                                  std::sin(2.0 * M_PI * x) *
                                  std::sin(2.0 * M_PI * y) *
                                  std::cos(M_PI * time / TT);
                });
        }
        u_mac.FillBoundary(geom[lev].periodicity());
        v_mac.FillBoundary(geom[lev].periodicity());
        w_mac.FillBoundary(geom[lev].periodicity());
    }
}

// Nothing to do afterward. Cell-centered velocity should be untouched for the
// sake of the test, even though this means the reported CFL will be
// meaningless.
void VortexPatchScalarVel::post_advance_work() {}

} // namespace amr_wind
