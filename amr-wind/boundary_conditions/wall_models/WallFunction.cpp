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
    amrex::Real mu;
    amrex::Real rho{1.0};
    amrex::Real amplitude;
    amrex::Real wavenumber;
    amrex::Real omega; 

    {
        amrex::ParmParse pp("BodyForce");
        amrex::Vector<amrex::Real> body_force{0.0, 0.0, 0.0};
        pp.getarr("magnitude", body_force);
        m_log_law.utau_mean = std::sqrt(std::sqrt(
            body_force[0] * body_force[0] + body_force[1] * body_force[1]));
    }
    {
        amrex::ParmParse pp("transport");
        pp.get("viscosity", mu);
    }
    {
        amrex::ParmParse pp("incflo");
        pp.query("density", rho);
    }
    m_log_law.nu = mu / rho;
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
    ///Dyanmic wave stress
    {
    amrex::ParmParse pp("wave_mosd"); // "wave_mosd" is the prefix used in the input file
    pp.query("amplitude", amplitude);
    pp.query("wavenumber", wavenumber);
    pp.query("frequency", omega);
    }
    m_mosd.amplitude = amplitude;
    m_mosd.wavenumber = wavenumber;
    m_mosd.omega = omega;


}

VelWallFunc::VelWallFunc(Field& /*unused*/, WallFunction& wall_func)
    : m_wall_func(wall_func)
{
    amrex::ParmParse pp("WallFunction");
    pp.query("wall_shear_stress_type", m_wall_shear_stress_type);
    m_wall_shear_stress_type = amrex::toLower(m_wall_shear_stress_type);

    if (m_wall_shear_stress_type == "constant" ||
        m_wall_shear_stress_type == "log_law" ||
        m_wall_shear_stress_type == "schumann" ||
	m_wall_shear_stress_type == "mosd")  {
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
        const auto& vold_lev = velocity.state(FieldState::Old)(lev);
        const auto& eta_lev = viscosity(lev);

	const auto& problo = geom.ProbLoArray();
        const auto& dx = geom.CellSizeArray();

        if (amrex::Gpu::notInLaunchRegion()) {
            mfi_info.SetDynamic(true);
        }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(vel_lev, mfi_info); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.validbox();
            const auto& varr = vel_lev.array(mfi);
            const auto& vold_arr = vold_lev.const_array(mfi);
            const auto& den = rho_lev.const_array(mfi);
            const auto& eta = eta_lev.const_array(mfi);

            ///Dyanmic wave stress
            // Getting the closest integer of the vertical grid point of value dx (z=dx)
            // this should be interpolated to be more accurate
            const amrex::Real dx_dz = dx[0] / dx[2];
            const int z_dx_integer = static_cast<int>(std::round(dx_dz));

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
                        //varr(i, j, k - 1, 2) = 0.0;
                        // Shear stress BC
                        //varr(i, j, k - 1, 0) =
                        //    tau.get_shear(uu, wspd) / mu * den(i, j, k);
                        //varr(i, j, k - 1, 1) =
                        //    tau.get_shear(vv, wspd) / mu * den(i, j, k);

                        ///Dyanmic wave stress
                        //if (m_wall_shear_stress_type == "mosd")
                        // {
                                                   
                        const amrex::Real xc = problo[0] + (i + 0.5) * dx[0];
			//const amrex::Real u_dx = vold_arr(i, j, k + z_dx_integer, 0);  
                        //const amrex::Real v_dx = vold_arr(i, j, k + z_dx_integer, 1);

                        const int z_dx_low = static_cast<int>(std::floor(dx[0] / dx[2]));
                        const int z_dx_up = z_dx_low + 1;
                        const amrex::Real z_diff = dx[0]  - (z_dx_low + 0.5) * dx[2]; // Difference between position and lower grid point

                        // Interpolate u and v using the lower and upper grid points
                        const amrex::Real u_low = vold_arr(i, j, z_dx_low, 0);
                        const amrex::Real u_up = vold_arr(i, j, z_dx_up, 0);
                        const amrex::Real v_low = vold_arr(i, j, z_dx_low, 1);
                        const amrex::Real v_up = vold_arr(i, j, z_dx_up, 1);

                       // Linear interpolation for u and v
                       const amrex::Real u_dx = u_low + (u_up - u_low) * (z_diff / dx[2]);
                       const amrex::Real v_dx = v_low + (v_up - v_low) * (z_diff / dx[2]);
                        varr(i, j, k - 1, 0) =
                            tau.get_shear(uu, wspd, u_dx, v_dx, xc, 0) / mu * den(i, j, k);
                        varr(i, j, k - 1, 1) =
                            tau.get_shear(vv, wspd, u_dx, v_dx, xc, 1) / mu * den(i, j, k);
                        varr(i, j, k - 1, 2) = 0.0;
                        
			//}
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
                        
			const amrex::Real xc = problo[0] + (i + 0.5) * dx[0];
                        const amrex::Real u_dx = vold_arr(i, j, k - 1 - z_dx_integer, 0);
			const amrex::Real v_dx = vold_arr(i, j, k - 1 - z_dx_integer, 1);
			
			// Dirichlet BC
                        varr(i, j, k, 2) = 0.0;

                        // Shear stress BC
                        varr(i, j, k, 0) =
                            -tau.get_shear(uu, wspd, u_dx, v_dx, xc, 0) / mu * den(i, j, k);
                        varr(i, j, k, 1) =
                            -tau.get_shear(vv, wspd, u_dx, v_dx, xc, 1) / mu * den(i, j, k);
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
            const auto& varr = vel_lev.array(mfi);
            const auto& den = rho_lev.const_array(mfi);
            const auto& eta = eta_lev.const_array(mfi);

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
    } else if (m_wall_shear_stress_type == "mosd") {
        m_wall_func.update_umean();
        m_wall_func.update_utau_mean();
        m_wall_func.update_time();
        auto tau = SimpleShearMOSD(m_wall_func.log_law(), m_wall_func.mosd());
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

void WallFunction::update_time()
{
    m_mosd.time = m_sim.time().current_time();

}
} // namespace amr_wind
