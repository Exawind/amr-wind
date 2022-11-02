#include "amr-wind/incflo.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"

#include <cmath>
#include <limits>

using namespace amrex;

/** Estimate the new timestep for adaptive timestepping algorithm
 *
 *  \param explicit_diffusion Flag indicating whether user has selected explicit
 *  treatment of diffusion term.
 *
 * Compute new \f$\Delta t\f$ by using the formula derived in "A Boundary
 * Condition Capturing Method for Multiphase Incompressible Flow" by Kang et al.
 * (JCP).
 *
 *  \f{align}
 *  \text{CFL} &= \frac{\Delta t}{2}  \left[\left(C+V \right) +
 * \sqrt{\left(C+V\right)^2 + 4 \left(\frac{|F_x|}{\Delta x} +
 * \frac{|F_y|}{\Delta y} + \frac{|F_z|}{\Delta
 * z} \right)} \right] && \\
 *  C &= \text{max} \left(\frac{|U_x|}{\Delta x} ,\ \frac{|U_y|}{\Delta y} ,\
 * \frac{|U_z|}{\Delta z}\right) && \text{Convection} \\ V &= 2
 * \left(\frac{\mu}{\rho}\right)_\mathrm{max} \left[\frac{1}{\Delta x^2} +
 * \frac{1}{\Delta y^2} + \frac{1}{\Delta z^2} \right] && \text{Diffusion} \f}
 *
 *  \f$F_x, F_y, F_z\f$ are acceleration due to forcing terms. By default,
 *  \f$C\f$ is always used for computing CFL. The term \f$V\f$ is only used when
 *  the user has chosen implicit or Crank-Nicholson scheme for diffusion, and
 *  contributions from forcing term when `time.use_force_cfl` is `true` (default
 *  is `true`).
 *
 */
void incflo::ComputeDt(bool explicit_diffusion)
{
    BL_PROFILE("amr-wind::incflo::ComputeDt");

    Real conv_cfl = 0.0;
    Real diff_cfl = 0.0;
    Real force_cfl = 0.0;
    const bool mesh_mapping = m_sim.has_mesh_mapping();

    const auto& den = density();
    amr_wind::Field const* mesh_fac =
        mesh_mapping
            ? &(m_repo.get_mesh_mapping_field(amr_wind::FieldLoc::CELL))
            : nullptr;

    for (int lev = 0; lev <= finest_level; ++lev) {
        auto const dxinv = geom[lev].InvCellSizeArray();
        MultiFab const& vel = icns().fields().field(lev);
        MultiFab const& vel_force = icns().fields().src_term(lev);
        MultiFab const& mu = icns().fields().mueff(lev);
        MultiFab const& rho = den(lev);

        auto const& vel_arr = vel.const_arrays();
        MultiArray4<Real const> fac_arr =
            mesh_mapping ? ((*mesh_fac)(lev).const_arrays())
                         : MultiArray4<Real const>();

        Real conv_lev = 0.0;
        Real mphase_conv_lev = 0.0;
        Real diff_lev = 0.0;
        Real force_lev = 0.0;

        conv_lev += amrex::ParReduce(
            TypeList<ReduceOpMax>{}, TypeList<Real>{}, vel, IntVect(0),
            [=] AMREX_GPU_HOST_DEVICE(
                int box_no, int i, int j, int k) -> GpuTuple<Real> {
                auto const& v_bx = vel_arr[box_no];

                amrex::Real fac_x =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 0)) : 1.0;
                amrex::Real fac_y =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 1)) : 1.0;
                amrex::Real fac_z =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 2)) : 1.0;

                return amrex::max(
                    amrex::Math::abs(v_bx(i, j, k, 0)) * dxinv[0] / fac_x,
                    amrex::Math::abs(v_bx(i, j, k, 1)) * dxinv[1] / fac_y,
                    amrex::Math::abs(v_bx(i, j, k, 2)) * dxinv[2] / fac_z,
                    -1.0);
            });

        if (m_sim.pde_manager().has_pde("VOF")) {
            MultiFab const& vof = m_repo.get_field("vof")(lev);
            auto const& vof_arr = vof.const_arrays();
            mphase_conv_lev += amrex::ParReduce(
                TypeList<ReduceOpMax>{}, TypeList<Real>{}, vel, IntVect(0),
                [=] AMREX_GPU_HOST_DEVICE(
                    int box_no, int i, int j, int k) -> GpuTuple<Real> {
                    auto const& v_bx = vel_arr[box_no];
                    auto const& vof_bx = vof_arr[box_no];

                    // Check for interface
                    auto is_near =
                        amr_wind::multiphase::interface_band(i, j, k, vof_bx);

                    // CFL calculation is not needed away from interface
                    amrex::Real result = 0.0;
                    if (is_near) {
                        // Near interface, evaluate CFL by sum of velocities
                        amrex::Real fac_x =
                            mesh_mapping ? (fac_arr[box_no](i, j, k, 0)) : 1.0;
                        amrex::Real fac_y =
                            mesh_mapping ? (fac_arr[box_no](i, j, k, 1)) : 1.0;
                        amrex::Real fac_z =
                            mesh_mapping ? (fac_arr[box_no](i, j, k, 2)) : 1.0;

                        result = amrex::Math::abs(v_bx(i, j, k, 0)) * dxinv[0] /
                                     fac_x +
                                 amrex::Math::abs(v_bx(i, j, k, 1)) * dxinv[1] /
                                     fac_y +
                                 amrex::Math::abs(v_bx(i, j, k, 2)) * dxinv[2] /
                                     fac_z;
                    }
                    return result;
                });
        }
        conv_lev = amrex::max(conv_lev, mphase_conv_lev);

        if (explicit_diffusion) {
            auto const& mu_arr = mu.const_arrays();
            auto const& rho_arr = rho.const_arrays();
            diff_lev += amrex::ParReduce(
                TypeList<ReduceOpMax>{}, TypeList<Real>{}, rho, IntVect(0),
                [=] AMREX_GPU_HOST_DEVICE(
                    int box_no, int i, int j, int k) -> GpuTuple<Real> {
                    auto const& mu_bx = mu_arr[box_no];
                    auto const& rho_bx = rho_arr[box_no];

                    amrex::Real fac_x =
                        mesh_mapping ? (fac_arr[box_no](i, j, k, 0)) : 1.0;
                    amrex::Real fac_y =
                        mesh_mapping ? (fac_arr[box_no](i, j, k, 1)) : 1.0;
                    amrex::Real fac_z =
                        mesh_mapping ? (fac_arr[box_no](i, j, k, 2)) : 1.0;

                    const Real dxinv2 =
                        2.0 * (dxinv[0] / fac_x * dxinv[0] / fac_x +
                               dxinv[1] / fac_y * dxinv[1] / fac_y +
                               dxinv[2] / fac_z * dxinv[2] / fac_z);

                    return amrex::max(
                        mu_bx(i, j, k) * dxinv2 / rho_bx(i, j, k), -1.0);
                });
        }

        if (m_time.use_force_cfl()) {
            auto const& vf_arr = vel_force.const_arrays();
            auto const& rho_arr = rho.const_arrays();
            force_lev += amrex::ParReduce(
                TypeList<ReduceOpMax>{}, TypeList<Real>{}, vel_force,
                IntVect(0),
                [=] AMREX_GPU_HOST_DEVICE(
                    int box_no, int i, int j, int k) -> GpuTuple<Real> {
                    auto const& vf_bx = vf_arr[box_no];
                    auto const& rho_bx = rho_arr[box_no];

                    amrex::Real fac_x =
                        mesh_mapping ? (fac_arr[box_no](i, j, k, 0)) : 1.0;
                    amrex::Real fac_y =
                        mesh_mapping ? (fac_arr[box_no](i, j, k, 1)) : 1.0;
                    amrex::Real fac_z =
                        mesh_mapping ? (fac_arr[box_no](i, j, k, 2)) : 1.0;

                    return amrex::max(
                        amrex::Math::abs(vf_bx(i, j, k, 0)) * dxinv[0] / fac_x /
                            rho_bx(i, j, k),
                        amrex::Math::abs(vf_bx(i, j, k, 1)) * dxinv[1] / fac_y /
                            rho_bx(i, j, k),
                        amrex::Math::abs(vf_bx(i, j, k, 2)) * dxinv[2] / fac_z /
                            rho_bx(i, j, k),
                        -1.0);
                });
        }

        conv_cfl = amrex::max(conv_cfl, conv_lev);
        diff_cfl = amrex::max(diff_cfl, diff_lev);
        force_cfl = amrex::max(force_cfl, force_lev);
    }

    ParallelAllReduce::Max<Real>(conv_cfl, ParallelContext::CommunicatorSub());
    if (explicit_diffusion) {
        ParallelAllReduce::Max<Real>(
            diff_cfl, ParallelContext::CommunicatorSub());
    }
    if (m_time.use_force_cfl()) {
        ParallelAllReduce::Max<Real>(
            force_cfl, ParallelContext::CommunicatorSub());
    }

    m_time.set_current_cfl(conv_cfl, diff_cfl, force_cfl);
}

void incflo::ComputePrescribeDt()
{
    BL_PROFILE("amr-wind::incflo::ComputePrescribeDt");

    Real conv_cfl = 0.0;
    const bool mesh_mapping = m_sim.has_mesh_mapping();

    amr_wind::Field const* mesh_fac =
        mesh_mapping
            ? &(m_repo.get_mesh_mapping_field(amr_wind::FieldLoc::CELL))
            : nullptr;

    for (int lev = 0; lev <= finest_level; ++lev) {
        auto const dxinv = geom[lev].InvCellSizeArray();
        auto const& uf_arr = m_repo.get_field("u_mac")(lev).const_arrays();
        auto const& vf_arr = m_repo.get_field("v_mac")(lev).const_arrays();
        auto const& wf_arr = m_repo.get_field("w_mac")(lev).const_arrays();

        MultiArray4<Real const> fac_arr =
            mesh_mapping ? ((*mesh_fac)(lev).const_arrays())
                         : MultiArray4<Real const>();

        Real conv_lev = 0.0;
        Real mphase_conv_lev = 0.0;

        conv_lev += amrex::ParReduce(
            TypeList<ReduceOpMax>{}, TypeList<Real>{},
            icns().fields().field(lev), IntVect(0),
            [=] AMREX_GPU_HOST_DEVICE(
                int box_no, int i, int j, int k) -> GpuTuple<Real> {
                auto const& umac = uf_arr[box_no];
                auto const& vmac = vf_arr[box_no];
                auto const& wmac = wf_arr[box_no];

                amrex::Real fac_x =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 0)) : 1.0;
                amrex::Real fac_y =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 1)) : 1.0;
                amrex::Real fac_z =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 2)) : 1.0;

                return amrex::max(
                    amrex::max(
                        amrex::Math::abs(umac(i, j, k)),
                        amrex::Math::abs(umac(i + 1, j, k))) *
                        dxinv[0] / fac_x,
                    amrex::max(
                        amrex::Math::abs(vmac(i, j, k)),
                        amrex::Math::abs(vmac(i, j + 1, k))) *
                        dxinv[1] / fac_y,
                    amrex::max(
                        amrex::Math::abs(wmac(i, j, k)),
                        amrex::Math::abs(wmac(i, j, k + 1))) *
                        dxinv[2] / fac_z,
                    -1.0);
            });

        if (m_sim.pde_manager().has_pde("VOF")) {
            MultiFab const& vof = m_repo.get_field("vof")(lev);
            auto const& vof_arr = vof.const_arrays();
            mphase_conv_lev += amrex::ParReduce(
                TypeList<ReduceOpMax>{}, TypeList<Real>{},
                icns().fields().field(lev), IntVect(0),
                [=] AMREX_GPU_HOST_DEVICE(
                    int box_no, int i, int j, int k) -> GpuTuple<Real> {
                    auto const& vof_bx = vof_arr[box_no];
                    auto const& umac = uf_arr[box_no];
                    auto const& vmac = vf_arr[box_no];
                    auto const& wmac = wf_arr[box_no];

                    // Check for interface
                    auto is_near =
                        amr_wind::multiphase::interface_band(i, j, k, vof_bx);

                    // CFL calculation is not needed away from interface
                    amrex::Real result = 0.0;
                    if (is_near) {
                        // Near interface, evaluate CFL by sum of velocities
                        amrex::Real fac_x =
                            mesh_mapping ? (fac_arr[box_no](i, j, k, 0)) : 1.0;
                        amrex::Real fac_y =
                            mesh_mapping ? (fac_arr[box_no](i, j, k, 1)) : 1.0;
                        amrex::Real fac_z =
                            mesh_mapping ? (fac_arr[box_no](i, j, k, 2)) : 1.0;

                        result = amrex::max(
                                     amrex::Math::abs(umac(i, j, k)),
                                     amrex::Math::abs(umac(i + 1, j, k))) *
                                     dxinv[0] / fac_x +
                                 amrex::max(
                                     amrex::Math::abs(vmac(i, j, k)),
                                     amrex::Math::abs(vmac(i, j + 1, k))) *
                                     dxinv[1] / fac_y +
                                 amrex::max(
                                     amrex::Math::abs(wmac(i, j, k)),
                                     amrex::Math::abs(wmac(i, j, k + 1))) *
                                     dxinv[2] / fac_z;
                    }
                    return result;
                });
        }
        conv_lev = amrex::max(conv_lev, mphase_conv_lev);

        conv_cfl = amrex::max(conv_cfl, conv_lev);
    }

    ParallelAllReduce::Max<Real>(conv_cfl, ParallelContext::CommunicatorSub());

    m_time.set_current_cfl(conv_cfl, 0.0, 0.0);
}
