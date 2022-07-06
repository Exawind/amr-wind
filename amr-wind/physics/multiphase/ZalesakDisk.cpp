#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/physics/multiphase/ZalesakDisk.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/core/field_ops.H"

namespace amr_wind {

ZalesakDisk::ZalesakDisk(CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
    , m_density(sim.repo().get_field("density"))
{
    amrex::ParmParse pp(identifier());
    pp.queryarr("location", m_loc, 0, AMREX_SPACEDIM);
    pp.query("radius", m_radius);
    pp.query("period", m_TT);
    amrex::ParmParse pinc("incflo");
    pinc.add("prescribe_velocity", true);
}

/** Initialize the velocity and levelset fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::ZalesakDiskFieldInit
 */
void ZalesakDisk::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& levelset = m_levelset(level);
    auto& density = m_density(level);
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();

    const auto& mphase = m_sim.physics_manager().get<MultiPhase>();
    const amrex::Real rho1 = mphase.rho1();
    const amrex::Real rho2 = mphase.rho2();

    auto& u_mac = m_sim.repo().get_field("u_mac")(level);
    auto& v_mac = m_sim.repo().get_field("v_mac")(level);
    auto& w_mac = m_sim.repo().get_field("w_mac")(level);

    const amrex::Real xc = m_loc[0];
    const amrex::Real yc = m_loc[1];
    const amrex::Real zc = m_loc[2];
    const amrex::Real radius = m_radius;
    const amrex::Real TT = m_TT;
    const amrex::Real width = m_width;
    const amrex::Real depth = m_depth;

    for (amrex::MFIter mfi(levelset); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.growntilebox();
        auto uf = u_mac.array(mfi);
        auto vf = v_mac.array(mfi);
        auto wf = w_mac.array(mfi);
        auto vel = velocity.array(mfi);
        auto phi = levelset.array(mfi);
        auto rho = density.array(mfi);
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                uf(i, j, k) = 2.0 * M_PI / TT * (0.5 - y);
                vf(i, j, k) = 2.0 * M_PI / TT * (x - 0.5);
                wf(i, j, k) = 0.0;

                vel(i, j, k, 0) = 2.0 * M_PI / TT * (0.5 - y);
                vel(i, j, k, 1) = 2.0 * M_PI / TT * (x - 0.5);
                vel(i, j, k, 2) = 0.0;
                // First define the sphere

                phi(i, j, k) =
                    radius - std::sqrt(
                                 (x - xc) * (x - xc) + (y - yc) * (y - yc) +
                                 (z - zc) * (z - zc));
                // then the slot
                amrex::Real eps = std::cbrt(dx[0] * dx[1] * dx[2]);
                if (y - yc <= radius && y - yc >= radius - depth &&
                    std::abs(x - xc) <= width &&
                    std::sqrt(
                        (x - xc) * (x - xc) + (y - yc) * (y - yc) +
                        (z - zc) * (z - zc)) < radius + eps) {
                    amrex::Real d1;
                    if (x > xc) {
                        d1 = std::abs(xc + width - x);
                    } else {
                        d1 = std::abs(xc - width - x);
                    }
                    const amrex::Real d2 = std::abs(y - (yc + radius - depth));
                    const amrex::Real min_dist = amrex::min(d1, d2);

                    phi(i, j, k) = -min_dist;
                }

                amrex::Real smooth_heaviside;
                eps = std::cbrt(2. * dx[0] * dx[1] * dx[2]);
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
            });
    }
}

void ZalesakDisk::pre_advance_work()
{

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

                    uf(i, j, k) = 2.0 * M_PI / TT * (0.5 - y);
                    vf(i, j, k) = 2.0 * M_PI / TT * (x - 0.5);
                    wf(i, j, k) = 0.0;
                });
        }
        u_mac.FillBoundary(geom[lev].periodicity());
        v_mac.FillBoundary(geom[lev].periodicity());
        w_mac.FillBoundary(geom[lev].periodicity());
    }
}

void ZalesakDisk::post_advance_work()
{

    const int nlevels = m_sim.repo().num_active_levels();
    const auto& geom = m_sim.mesh().Geom();

    // Overriding the velocity field
    for (int lev = 0; lev < nlevels; ++lev) {
        for (amrex::MFIter mfi(m_velocity(lev)); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            const auto& dx = geom[lev].CellSizeArray();
            const auto& problo = geom[lev].ProbLoArray();
            const amrex::Real TT = m_TT;
            auto vel = m_velocity(lev).array(mfi);
            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];

                    vel(i, j, k, 0) = 2.0 * M_PI / TT * (0.5 - y);
                    vel(i, j, k, 1) = 2.0 * M_PI / TT * (x - 0.5);
                    vel(i, j, k, 2) = 0.0;
                });
        }
    }

    m_velocity.fillpatch(m_sim.time().current_time());
}

} // namespace amr_wind
