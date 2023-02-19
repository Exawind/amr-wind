#include "amr-wind/physics/VortexRing.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/trig_ops.H"
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_FillPatchUtil.H>
#include "amr-wind/core/FieldRepo.H"

namespace amr_wind {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real FatCore::operator()(
    const amrex::Real r,
    const amrex::Real /*unused*/,
    const amrex::Real z,
    const amrex::Real R,
    const amrex::Real Gamma,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const int /*unused*/,
    const int* /*unused*/,
    const amrex::Real* /*unused*/,
    const amrex::Real* /*unused*/) const
{
    amrex::Real Rsq = std::pow(R, 2);
    const amrex::Real ssq = std::pow(z, 2) + std::pow(r - R, 2);
    if (ssq < Rsq) {
        return 0.54857674 * Gamma / Rsq * std::exp(-4 * ssq / (Rsq - ssq));
    }
    return 0.0;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real CollidingRings::operator()(
    const amrex::Real r,
    const amrex::Real theta,
    const amrex::Real z,
    const amrex::Real R,
    const amrex::Real Gamma,
    const amrex::Real delta,
    const amrex::Real dz,
    const amrex::Real perturbation_amplitude,
    const int num_modes,
    const int* perturbation_modes,
    const amrex::Real* perturbation_phases_1,
    const amrex::Real* perturbation_phases_2) const
{
    amrex::Real dr1 = 0.0;
    for (int i = 0; i < num_modes; ++i) {
        dr1 +=
            perturbation_amplitude *
            std::cos(perturbation_modes[i] * theta - perturbation_phases_1[i]);
    }
    amrex::Real dr2 = 0.0;
    for (int i = 0; i < num_modes; ++i) {
        dr2 +=
            perturbation_amplitude *
            std::cos(perturbation_modes[i] * theta - perturbation_phases_2[i]);
    }
    amrex::Real vortheta_1 =
        -Gamma / (utils::pi() * std::pow(delta, 2)) *
        std::exp(
            -(std::pow(z + dz / 2, 2) + std::pow((r * (1 + dr1) - R), 2)) /
            std::pow(delta, 2));
    amrex::Real vortheta_2 =
        Gamma / (utils::pi() * std::pow(delta, 2)) *
        std::exp(
            -(std::pow(z - dz / 2, 2) + std::pow((r * (1 + dr2) - R), 2)) /
            std::pow(delta, 2));
    return vortheta_1 + vortheta_2;
}

VortexRing::VortexRing(const CFDSim& sim)
    : m_sim(sim)
    , m_repo(sim.repo())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{
    {
        amrex::ParmParse pp("incflo");
        pp.query("density", m_rho);
    }

    {
        amrex::ParmParse pp("vortexring");
        pp.query("type", m_vortexringtype);
        pp.query("R", m_R);
        pp.query("Gamma", m_Gamma);
        pp.query("delta", m_delta);
        pp.query("dz", m_dz);
        pp.query("perturbation_amplitude", m_perturbation_amplitude);
        pp.queryarr("perturbation_modes", m_perturbation_modes);
        pp.queryarr("perturbation_phases_1", m_perturbation_phases_1);
        pp.queryarr("perturbation_phases_2", m_perturbation_phases_2);
    }
}

/** Initialize the velocity and density fields at the beginning of the
 *  simulation.
 */
void VortexRing::initialize_fields(int level, const amrex::Geometry& /*geom*/)
{
    m_density(level).setVal(m_rho);
    m_velocity(level).setVal(0.0, 0, AMREX_SPACEDIM);
}

void VortexRing::post_init_actions()
{

    // only call for startup and not for restart
    if (m_sim.time().time_index() > 0) {
        return;
    }

    if (m_vortexringtype == "fatcore") {
        initialize_velocity(FatCore());
    } else if (m_vortexringtype == "collidingrings") {
        initialize_velocity(CollidingRings());
    }
}

template <typename VortexRingType>
void VortexRing::initialize_velocity(const VortexRingType& vorticity_theta)
{

    const int nghost = 1;
    auto minusvorticity =
        m_repo.create_scratch_field(AMREX_SPACEDIM, nghost, FieldLoc::NODE);
    auto vectorpotential =
        m_repo.create_scratch_field(AMREX_SPACEDIM, nghost, FieldLoc::NODE);

    const amrex::Real R = m_R;
    const amrex::Real Gamma = m_Gamma;
    const amrex::Real delta = m_delta;
    const amrex::Real dz = m_dz;
    const amrex::Real perturbation_amplitude = m_perturbation_amplitude;

    amrex::Gpu::DeviceVector<int> perturbation_modes_d;
    amrex::Gpu::DeviceVector<amrex::Real> perturbation_phases_1_d;
    amrex::Gpu::DeviceVector<amrex::Real> perturbation_phases_2_d;
    const int num_modes = static_cast<int>(m_perturbation_modes.size());
    if (num_modes > 0) {

        perturbation_modes_d.resize(num_modes);
        perturbation_phases_1_d.resize(num_modes);
        perturbation_phases_2_d.resize(num_modes);

        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_perturbation_modes.begin(),
            m_perturbation_modes.end(), perturbation_modes_d.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_perturbation_phases_1.begin(),
            m_perturbation_phases_1.end(), perturbation_phases_1_d.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_perturbation_phases_2.begin(),
            m_perturbation_phases_2.end(), perturbation_phases_2_d.begin());
    }

    const int* pm = perturbation_modes_d.data();
    const amrex::Real* pp1 = perturbation_phases_1_d.data();
    const amrex::Real* pp2 = perturbation_phases_2_d.data();

    for (int level = 0; level <= m_repo.mesh().finestLevel(); ++level) {

        const auto& problo = m_repo.mesh().Geom(level).ProbLoArray();

        for (amrex::MFIter mfi(m_velocity(level)); mfi.isValid(); ++mfi) {
            const auto& dx = m_repo.mesh().Geom(level).CellSizeArray();
            const auto& nbx = mfi.nodaltilebox();
            auto minusvort = (*minusvorticity)(level).array(mfi);

            amrex::ParallelFor(
                nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + i * dx[0];
                    const amrex::Real y = problo[1] + j * dx[1];
                    const amrex::Real z = problo[2] + k * dx[2];
                    const amrex::Real r =
                        std::sqrt(std::pow(x, 2) + std::pow(y, 2));
                    const amrex::Real theta = std::atan2(y, x);
                    const amrex::Real vortheta = vorticity_theta(
                        r, theta, z, R, Gamma, delta, dz,
                        perturbation_amplitude, num_modes, pm, pp1, pp2);
                    minusvort(i, j, k, 0) = std::sin(theta) * vortheta;
                    minusvort(i, j, k, 1) = -std::cos(theta) * vortheta;
                    minusvort(i, j, k, 2) = 0.0;
                });
        }
    }

    amrex::LPInfo info;
    const auto& mesh = m_velocity.repo().mesh();

    amrex::MLNodeLaplacian linop(
        mesh.Geom(0, mesh.finestLevel()), mesh.boxArray(0, mesh.finestLevel()),
        mesh.DistributionMap(0, mesh.finestLevel()), info, {}, 1.0);

    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> bclo;
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> bchi;

    const auto& bctype = m_velocity.bc_type();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (mesh.Geom(0).isPeriodic(dir)) {
            bclo[dir] = amrex::LinOpBCType::Periodic;
            bchi[dir] = amrex::LinOpBCType::Periodic;
        } else {

            switch (bctype[amrex::Orientation(dir, amrex::Orientation::low)]) {
            case BC::pressure_inflow:
            case BC::pressure_outflow: {
                bclo[dir] = amrex::LinOpBCType::Dirichlet;
                break;
            }
            default:
                bclo[dir] = amrex::LinOpBCType::Neumann;
                break;
            };

            switch (bctype[amrex::Orientation(dir, amrex::Orientation::high)]) {
            case BC::pressure_inflow:
            case BC::pressure_outflow: {
                bchi[dir] = amrex::LinOpBCType::Dirichlet;
                break;
            }
            default:
                bchi[dir] = amrex::LinOpBCType::Neumann;
                break;
            };
        }
    }

    linop.setDomainBC(bclo, bchi);
    amrex::MLMG mlmg(linop);

    mlmg.setVerbose(2);
    const amrex::Real rel_tol = 1.0e-13;
    const amrex::Real abs_tol = 1.0e-13;

    for (int level = 0; level <= m_repo.mesh().finestLevel(); ++level) {
        (*vectorpotential)(level).setVal(0.0, 0, AMREX_SPACEDIM, 1);
    }

    // might be able to skip z-dir since vorticity is 0.0
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        auto vectorpot = (*vectorpotential).subview(i, 1);
        auto minusvort = (*minusvorticity).subview(i, 1);
        mlmg.solve(
            vectorpot.vec_ptrs(), minusvort.vec_const_ptrs(), rel_tol, abs_tol);
    }

    for (int level = 0; level <= m_repo.mesh().finestLevel(); ++level) {
        auto& velocity = m_velocity(level);

        for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();

            const auto& dxinv = m_repo.mesh().Geom(level).InvCellSizeArray();
            auto vel = velocity.array(mfi);
            auto psi = (*vectorpotential)(level).array(mfi);
            const amrex::Real facx = amrex::Real(0.25) * dxinv[0];
            const amrex::Real facy = amrex::Real(0.25) * dxinv[1];
            const amrex::Real facz = amrex::Real(0.25) * dxinv[2];

            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real dpsix_dy =
                        facy *
                        (-psi(i, j, k, 0) - psi(i + 1, j, k, 0) +
                         psi(i, j + 1, k, 0) + psi(i + 1, j + 1, k, 0) -
                         psi(i, j, k + 1, 0) - psi(i + 1, j, k + 1, 0) +
                         psi(i, j + 1, k + 1, 0) + psi(i + 1, j + 1, k + 1, 0));
                    const amrex::Real dpsix_dz =
                        facz *
                        (-psi(i, j, k, 0) - psi(i + 1, j, k, 0) -
                         psi(i, j + 1, k, 0) - psi(i + 1, j + 1, k, 0) +
                         psi(i, j, k + 1, 0) + psi(i + 1, j, k + 1, 0) +
                         psi(i, j + 1, k + 1, 0) + psi(i + 1, j + 1, k + 1, 0));
                    const amrex::Real dpsiy_dx =
                        facx *
                        (-psi(i, j, k, 1) + psi(i + 1, j, k, 1) -
                         psi(i, j + 1, k, 1) + psi(i + 1, j + 1, k, 1) -
                         psi(i, j, k + 1, 1) + psi(i + 1, j, k + 1, 1) -
                         psi(i, j + 1, k + 1, 1) + psi(i + 1, j + 1, k + 1, 1));
                    const amrex::Real dpsiy_dz =
                        facz *
                        (-psi(i, j, k, 1) - psi(i + 1, j, k, 1) -
                         psi(i, j + 1, k, 1) - psi(i + 1, j + 1, k, 1) +
                         psi(i, j, k + 1, 1) + psi(i + 1, j, k + 1, 1) +
                         psi(i, j + 1, k + 1, 1) + psi(i + 1, j + 1, k + 1, 1));
                    const amrex::Real dpsiz_dx =
                        facx *
                        (-psi(i, j, k, 2) + psi(i + 1, j, k, 2) -
                         psi(i, j + 1, k, 2) + psi(i + 1, j + 1, k, 2) -
                         psi(i, j, k + 1, 2) + psi(i + 1, j, k + 1, 2) -
                         psi(i, j + 1, k + 1, 2) + psi(i + 1, j + 1, k + 1, 2));
                    const amrex::Real dpsiz_dy =
                        facy *
                        (-psi(i, j, k, 2) - psi(i + 1, j, k, 2) +
                         psi(i, j + 1, k, 2) + psi(i + 1, j + 1, k, 2) -
                         psi(i, j, k + 1, 2) - psi(i + 1, j, k + 1, 2) +
                         psi(i, j + 1, k + 1, 2) + psi(i + 1, j + 1, k + 1, 2));

                    vel(i, j, k, 0) = -(dpsiz_dy - dpsiy_dz);
                    vel(i, j, k, 1) = -(dpsix_dz - dpsiz_dx);
                    vel(i, j, k, 2) = -(dpsiy_dx - dpsix_dy);
                });
        }
    }
}

} // namespace amr_wind
