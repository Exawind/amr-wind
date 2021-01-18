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
{
    amrex::ParmParse pp(identifier());
    pp.queryarr("location", m_loc, 0, AMREX_SPACEDIM);
    pp.query("radius", m_radius);
    pp.query("period", m_TT);
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
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();

    const amrex::Real xc = m_loc[0];
    const amrex::Real yc = m_loc[1];
    const amrex::Real zc = m_loc[2];
    const amrex::Real radius = m_radius;
    const amrex::Real TT = m_TT;
    const amrex::Real width = m_width;
    const amrex::Real depth = m_depth;

    for (amrex::MFIter mfi(levelset); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.growntilebox();
        auto vel = velocity.array(mfi);
        auto phi = levelset.array(mfi);
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                vel(i, j, k, 0) = 2.0 * M_PI / TT * (0.5 - y);
                vel(i, j, k, 1) = 2.0 * M_PI / TT * (x - 0.5);
                vel(i, j, k, 2) = 0.0;

                phi(i, j, k) =
                    radius - std::sqrt(
                                 (x - xc) * (x - xc) + (y - yc) * (y - yc) +
                                 (z - zc) * (z - zc));

                if (y - yc <= radius && y - yc > radius - depth &&
                    std::abs(x - xc) <= width &&
                    std::sqrt(
                        (x - xc) * (x - xc) + (y - yc) * (y - yc) +
                        (z - zc) * (z - zc)) <= radius) {
                    amrex::Real d1, d2;
                    if (x > xc) {
                        d1 = std::abs(xc + width - x);
                    } else {
                        d1 = std::abs(xc - width - x);
                    }

                    d2 = std::abs(y - (yc + radius - depth));
                    amrex::Real min_dist = amrex::min(d1, d2);

                    phi(i, j, k) = -min_dist;
                }
            });
    }
}

void ZalesakDisk::pre_advance_work() {}

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
