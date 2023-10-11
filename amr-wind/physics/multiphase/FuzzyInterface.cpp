#include "amr-wind/physics/multiphase/FuzzyInterface.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/core/field_ops.H"

namespace amr_wind {

FuzzyInterface::FuzzyInterface(CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_vof(sim.repo().get_field("vof"))
    , m_pressure(sim.repo().get_field("p"))
{
    amrex::ParmParse pp(identifier());
    pp.query("water_level", m_waterlevel);
    pp.get("lo_x_vofjump", m_lx_vofjump);
    pp.get("hi_x_vofjump", m_hx_vofjump);
    pp.query("interface_thickness", m_intf_th);

    amrex::ParmParse pp_multiphase("MultiPhase");
    pp_multiphase.add("water_level", m_waterlevel);

    auto& mphase = sim.physics_manager().get<MultiPhase>();
    m_rho1 = mphase.rho1();
    m_rho2 = mphase.rho2();
    m_gravity = mphase.gravity();

    // This physics case specifies vof directly
    mphase.set_intf_init_method_to_vof();
}

/** Initialize the velocity and levelset fields at the beginning of the
 *  simulation.
 *
 */
void FuzzyInterface::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    velocity.setVal(0.0);

    auto& volfrac = m_vof(level);
    auto& pressure = m_pressure(level);
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    const auto& probhi = geom.ProbHiArray();
    const amrex::Real water_level = m_waterlevel;
    const amrex::Real i_th = m_intf_th;
    const amrex::Real lx_vj = m_lx_vofjump;
    const amrex::Real hx_vj = m_hx_vofjump;
    const amrex::Real rho1 = m_rho1;
    const amrex::Real rho2 = m_rho2;
    const amrex::Real grav_z = m_gravity[2];

    for (amrex::MFIter mfi(volfrac); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        auto vof = volfrac.array(mfi);
        auto p = pressure.array(mfi);

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                const amrex::Real zbtm = problo[2] + (k + 0.0) * dx[2];
                const amrex::Real vof_sharp =
                    std::max(0.0, std::min(1.0, (water_level - zbtm) / dx[2]));
                const amrex::Real vof_smooth =
                    -0.5 * (std::erf((z - water_level) / m_intf_th) + 1.0) +
                    1.0;
                if (x < lx_vj || x > hx_vj) {
                    // Sharp vof
                    vof(i, j, k) = vof_sharp;
                } else {
                    vof(i, j, k) = vof_smooth;
                }
            });

        amrex::Box const& nbx = mfi.grownnodaltilebox();
        amrex::ParallelFor(
            nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // For pressure nodes, no offset
                const amrex::Real x = problo[0] + i * dx[0];
                const amrex::Real y = problo[1] + j * dx[1];
                const amrex::Real z = problo[2] + k * dx[2];

                // Sharp interpretation
                amrex::Real ih_g = amrex::max(
                    0.0, amrex::min(probhi[2] - water_level, probhi[2] - z));
                amrex::Real ih_l = amrex::max(
                    0.0, amrex::min(water_level - z, water_level - problo[2]));
                const amrex::Real irho = rho1 * ih_l + rho2 * ih_g;
                const amrex::Real psharp = -irho * grav_z;

                // Smooth interpretation
                const double x0 = (z - water_level) / i_th;
                const double x1 = (probhi[2] - water_level) / i_th;
                const double int_erf0 =
                    (x0 * std::erf(x0) +
                     std::exp(-std::pow(x0, 2)) / std::sqrt(M_PI)) *
                    i_th;
                const double int_erf1 =
                    (x1 * std::erf(x1) +
                     std::exp(-std::pow(x1, 2)) / std::sqrt(M_PI)) *
                    i_th;
                const double int_vof0 = 0.5 * z - 0.5 * int_erf0;
                const double int_vof1 = 0.5 * probhi[2] - 0.5 * int_erf1;
                const double int_rho0 = rho1 * int_vof0 + rho2 * (z - int_vof0);
                const double int_rho1 =
                    rho1 * int_vof1 + rho2 * (probhi[2] - int_vof1);

                // g * integral(rho)dz
                const amrex::Real psmooth = -(int_rho1 - int_rho0) * grav_z;

                // Populate pressure
                p(i, j, k) = (x < lx_vj || x > hx_vj) ? psharp : psmooth;
            });
    }
}

} // namespace amr_wind
