#include "amr-wind/boundary_conditions/wall_models/WallFunction.H"
#include "amr-wind/boundary_conditions/wall_models/ShearStressSimple.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/diffusion/diffusion.H"

#include <cmath>

#include "AMReX_ParmParse.H"
#include "AMReX_Print.H"
#include "AMReX_ParallelDescriptor.H"

namespace amr_wind {

WallFunction::WallFunction(CFDSim& sim)
    : m_sim(sim), m_mesh(m_sim.mesh()), m_pa_vel(sim, m_direction)
{
    {
        amrex::ParmParse pp("BodyForce");
        amrex::Vector<amrex::Real> body_force{{0.0, 0.0, 0.0}};
        pp.getarr("magnitude", body_force);
        m_log_law.utau_mean = std::sqrt(std::sqrt(
            body_force[0] * body_force[0] + body_force[1] * body_force[1]));
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::abs(body_force[2]) < 1e-16,
            "body force in z should be zero for this wall function");
    }
    {
        amrex::ParmParse pp("transport");
        pp.query("viscosity", m_log_law.nu);
    }
    {
        amrex::ParmParse pp("WallFunction");
        // Reference height for log-law
        if (pp.contains("log_law_ref_index")) {
            pp.get("log_law_ref_index", m_log_law.ref_index);
        }
        const auto& geom = m_mesh.Geom(0);
        m_log_law.zref =
            (geom.ProbLo(m_direction) +
             (m_log_law.ref_index + 0.5) * geom.CellSize(m_direction));
    }
}

VelWallFunc::VelWallFunc(Field& /*unused*/, WallFunction& wall_func)
    : m_wall_func(wall_func)
{
    amrex::ParmParse pp("WallFunction");
    pp.query("wall_shear_stress_type", m_wall_shear_stress_type);
    m_wall_shear_stress_type = amrex::toLower(m_wall_shear_stress_type);

    if (m_wall_shear_stress_type == "constant" ||
        m_wall_shear_stress_type == "log_law" ||
        m_wall_shear_stress_type == "schumann") {
        amrex::Print() << "Shear Stress model: " << m_wall_shear_stress_type
                       << std::endl;
    } else {
        amrex::Abort("Shear Stress wall model input mistake");
    }
}

template <typename ShearStressSimple>
void VelWallFunc::wall_model(
    Field& velocity, const FieldState rho_state, const ShearStressSimple& tau)
{
    BL_PROFILE("amr-wind::VelWallFunc");
    constexpr int idim = 2;
    const auto& repo = velocity.repo();
    const auto& density = repo.get_field("density", rho_state);
    const auto& viscosity = repo.get_field("velocity_mueff");
    const int nlevels = repo.num_active_levels();
    const auto idx_offset = tau.m_ll.ref_index;

    amrex::Orientation zlo(amrex::Direction::z, amrex::Orientation::low);
    amrex::Orientation zhi(amrex::Direction::z, amrex::Orientation::high);
    if ((velocity.bc_type()[zlo] != BC::wall_model) &&
        (velocity.bc_type()[zhi] != BC::wall_model)) {
        return;
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = repo.mesh().Geom(lev);
        const auto& domain = geom.Domain();
        amrex::MFItInfo mfi_info{};

        const auto& rho_lev = density(lev);
        auto& vel_lev = velocity(lev);
        auto& vold_lev = velocity.state(FieldState::Old)(lev);
        const auto& eta_lev = viscosity(lev);

        if (amrex::Gpu::notInLaunchRegion()) {
            mfi_info.SetDynamic(true);
        }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(vel_lev, mfi_info); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.validbox();
            auto varr = vel_lev.array(mfi);
            auto vold_arr = vold_lev.array(mfi);
            auto den = rho_lev.array(mfi);
            auto eta = eta_lev.array(mfi);

            if (bx.smallEnd(idim) == domain.smallEnd(idim) &&
                velocity.bc_type()[zlo] == BC::wall_model) {
                amrex::ParallelFor(
                    amrex::bdryLo(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real mu = eta(i, j, k);
                        const amrex::Real uu =
                            vold_arr(i, j, k + idx_offset, 0);
                        const amrex::Real vv =
                            vold_arr(i, j, k + idx_offset, 1);
                        const amrex::Real wspd = std::sqrt(uu * uu + vv * vv);
                        // Dirichlet BC
                        varr(i, j, k - 1, 2) = 0.0;

                        // Shear stress BC
                        varr(i, j, k - 1, 0) =
                            tau.get_shear(uu, wspd) / mu * den(i, j, k);
                        varr(i, j, k - 1, 1) =
                            tau.get_shear(vv, wspd) / mu * den(i, j, k);
                    });
            }

            if (bx.bigEnd(idim) == domain.bigEnd(idim) &&
                velocity.bc_type()[zhi] == BC::wall_model) {
                amrex::ParallelFor(
                    amrex::bdryHi(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real mu = eta(i, j, k - 1);
                        const amrex::Real uu =
                            vold_arr(i, j, k - 1 - idx_offset, 0);
                        const amrex::Real vv =
                            vold_arr(i, j, k - 1 - idx_offset, 1);
                        const amrex::Real wspd = std::sqrt(uu * uu + vv * vv);
                        // Dirichlet BC
                        varr(i, j, k, 2) = 0.0;

                        // Shear stress BC
                        varr(i, j, k, 0) =
                            -tau.get_shear(uu, wspd) / mu * den(i, j, k);
                        varr(i, j, k, 1) =
                            -tau.get_shear(vv, wspd) / mu * den(i, j, k);
                    });
            }
        }
    }
}

void VelWallFunc::wall_model(
    Field& velocity, const FieldState rho_state, const amrex::Real utau)
{
    BL_PROFILE("amr-wind::VelWallFunc");
    constexpr int idim = 2;
    const auto& repo = velocity.repo();
    const auto& density = repo.get_field("density", rho_state);
    const auto& viscosity = repo.get_field("velocity_mueff");
    const int nlevels = repo.num_active_levels();

    amrex::Orientation zlo(amrex::Direction::z, amrex::Orientation::low);
    amrex::Orientation zhi(amrex::Direction::z, amrex::Orientation::high);
    if ((velocity.bc_type()[zlo] != BC::wall_model) &&
        (velocity.bc_type()[zhi] != BC::wall_model)) {
        return;
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = repo.mesh().Geom(lev);
        const auto& domain = geom.Domain();
        amrex::MFItInfo mfi_info{};

        const auto& rho_lev = density(lev);
        auto& vel_lev = velocity(lev);
        const auto& eta_lev = viscosity(lev);

        if (amrex::Gpu::notInLaunchRegion()) {
            mfi_info.SetDynamic(true);
        }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(vel_lev, mfi_info); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.validbox();
            auto varr = vel_lev.array(mfi);
            auto den = rho_lev.array(mfi);
            auto eta = eta_lev.array(mfi);

            if (bx.smallEnd(idim) == domain.smallEnd(idim) &&
                velocity.bc_type()[zlo] == BC::wall_model) {
                amrex::ParallelFor(
                    amrex::bdryLo(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real mu = eta(i, j, k);

                        // Dirichlet BC
                        varr(i, j, k - 1, 2) = 0.0;

                        // Shear stress BC
                        varr(i, j, k - 1, 0) = utau * utau / mu * den(i, j, k);
                        varr(i, j, k - 1, 1) = 0.0;
                    });
            }

            if (bx.bigEnd(idim) == domain.bigEnd(idim) &&
                velocity.bc_type()[zhi] == BC::wall_model) {
                amrex::ParallelFor(
                    amrex::bdryHi(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real mu = eta(i, j, k - 1);

                        // Dirichlet BC
                        varr(i, j, k, 2) = 0.0;

                        // Shear stress BC
                        varr(i, j, k, 0) = utau * utau / mu * den(i, j, k);
                        varr(i, j, k, 1) = 0.0;
                    });
            }
        }
    }
}

void VelWallFunc::operator()(Field& velocity, const FieldState rho_state)
{
    if (m_wall_shear_stress_type == "constant") {
        wall_model(velocity, rho_state, m_wall_func.utau());
    } else if (m_wall_shear_stress_type == "log_law") {
        m_wall_func.update_umean();
        m_wall_func.update_utau_mean();
        auto tau = SimpleShearLogLaw(m_wall_func.log_law());
        wall_model(velocity, rho_state, tau);
    } else if (m_wall_shear_stress_type == "schumann") {
        m_wall_func.update_umean();
        auto tau = SimpleShearSchumann(m_wall_func.log_law());
        wall_model(velocity, rho_state, tau);
    }
}

void WallFunction::update_umean()
{
    m_pa_vel();
    m_log_law.wspd_mean =
        m_pa_vel.line_hvelmag_average_interpolated(m_log_law.zref);
}

void WallFunction::update_utau_mean() { m_log_law.update_utau_mean(); }
} // namespace amr_wind
