#include "amr-wind/physics/ChannelFlow.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/DirectionSelector.H"

namespace amr_wind {
namespace channel_flow {

ChannelFlow::ChannelFlow(CFDSim& sim)
    : m_time(sim.time()), m_repo(sim.repo()), m_mesh(sim.mesh())
{

    {
        amrex::ParmParse pp("ChannelFlow");
        pp.query("normal_direction", m_norm_dir);

        pp.query("density", m_rho);
        pp.query("re_tau", m_re_tau);
        pp.query("tke0", m_tke0);
        pp.query("sdr0", m_sdr0);
    }
    {
        amrex::Real mu;
        amrex::ParmParse pp("transport");
        pp.query("viscosity", mu);
        // Assumes a boundary layer height of 1.0
        m_utau = mu * m_re_tau / (m_rho * 1.0);
        m_ytau = mu / (m_utau * m_rho);
    }
}

/** Initialize the velocity, density, tke and sdr fields at the beginning of the
 *  simulation.
 */
void ChannelFlow::initialize_fields(int level, const amrex::Geometry& geom)
{

    switch (m_norm_dir) {
    case 1:
        initialize_fields(level, geom, YDir(), 1);
        break;
    case 2:
        initialize_fields(level, geom, ZDir(), 2);
        break;
    default:
        amrex::Abort("axis must be equal to 1 or 2");
        break;
    }
}

template <typename IndexSelector>
void ChannelFlow::initialize_fields(
    int level,
    const amrex::Geometry& geom,
    const IndexSelector& idxOp,
    const int n_idx)
{

    const amrex::Real kappa = m_kappa;
    const amrex::Real y_tau = m_ytau;
    const amrex::Real utau = m_utau;
    auto& velocity = m_repo.get_field("velocity")(level);
    auto& density = m_repo.get_field("density")(level);
    auto& tke = m_repo.get_field("tke")(level);
    auto& sdr = m_repo.get_field("sdr")(level);
    auto& walldist = m_repo.get_field("wall_dist")(level);

    density.setVal(m_rho);
    tke.setVal(m_tke0);
    sdr.setVal(m_sdr0);

    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        const auto& dx = geom.CellSizeArray();
        const auto& problo = geom.ProbLoArray();
        auto vel = velocity.array(mfi);
        auto wd = walldist.array(mfi);

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const int n_ind = idxOp(i, j, k);
                amrex::Real h = problo[n_idx] + (n_ind + 0.5) * dx[n_idx];
                if (h > 1.0) h = 2.0 - h;
                wd(i, j, k) = h;
                const amrex::Real hp = h / y_tau;
                vel(i, j, k, 0) =
                    utau * (1. / kappa * std::log1p(kappa * hp) +
                            7.8 * (1.0 - std::exp(-hp / 11.0) -
                                   (hp / 11.0) * std::exp(-hp / 3.0)));
                // vel(i,j,k,0) = 22.0;
                vel(i, j, k, 1) = 0.0;
                vel(i, j, k, 2) = 0.0;
            });
    }
}

void ChannelFlow::post_init_actions() {}

void ChannelFlow::post_advance_work() {}

} // namespace channel_flow
} // namespace amr_wind
