#include <incflo.H>
#include "Physics.H"

using namespace amrex;

void incflo::compute_vel_forces (Vector<MultiFab*> const& vel_forces,
                                 Vector<MultiFab const*> const& velocity,
                                 Vector<MultiFab const*> const& density,
                                 Vector<MultiFab const*> const& tracer)
{
    // FIXME: Clean up problem type specific logic
    if (m_probtype != 5) {
        for (int lev=0; lev <= finest_level; ++lev) {
            compute_vel_pressure_terms(lev, *vel_forces[lev], *density[lev]);

            for (auto& pp: m_physics) {
                pp->add_momentum_sources(
                    Geom(lev), *density[lev], *velocity[lev], *tracer[lev], *vel_forces[lev]);
            }
        }
    } else {
        // fixme this is the non boussinesq pathway and this could go away except that rayleigh-taylor needs it
        // there might be a bug in here since gravity is added twice
        // it's added twice because background pressure used to put it in m_gp0
        // need verification to know for sure
        // vel_f += - dpdx/rho + g*rho_0/rho + g
        for (int lev = 0; lev <= finest_level; ++lev)
            compute_vel_forces_on_level(
                lev, *vel_forces[lev], *density[lev], *tracer[lev]);
    }
}

void incflo::compute_vel_pressure_terms(int lev, amrex::MultiFab& vel_forces,
                                        const amrex::MultiFab& density)
{
    const amrex::GpuArray<amrex::Real, 3> l_gp0{m_gp0[0], m_gp0[1], m_gp0[2]};

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(vel_forces, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.tilebox();
        Array4<Real> const& vel_f = vel_forces.array(mfi);
        Array4<Real const> const& rho = density.const_array(mfi);
        Array4<Real const> const& gradp = grad_p()(lev).const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            amrex::Real rhoinv = 1.0 / rho(i, j, k);
            vel_f(i, j, k, 0) = -(gradp(i, j, k, 0) + l_gp0[0]) * rhoinv;
            vel_f(i, j, k, 1) = -(gradp(i, j, k, 1) + l_gp0[1]) * rhoinv;
            vel_f(i, j, k, 2) = -(gradp(i, j, k, 2) + l_gp0[2]) * rhoinv;
        });
    }
}

void incflo::compute_vel_forces_on_level (int lev,
                                                MultiFab& vel_forces,
                                          const MultiFab& density,
                                          const MultiFab& /* tracer */)
{
    GpuArray<Real,3> l_gravity{m_gravity[0],m_gravity[1],m_gravity[2]};
    GpuArray<Real,3> l_gp0{m_gp0[0] - m_gravity[0] * m_ro_0,
                           m_gp0[1] - m_gravity[1] * m_ro_0,
                           m_gp0[2] - m_gravity[2] * m_ro_0};

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(vel_forces,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const& vel_f = vel_forces.array(mfi);
        Array4<Real const> const& rho = density.const_array(mfi);
        Array4<Real const> const& gradp = grad_p()(lev).const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoinv = 1.0/rho(i,j,k);
            vel_f(i,j,k,0) = -(gradp(i,j,k,0)+l_gp0[0])*rhoinv + l_gravity[0];
            vel_f(i,j,k,1) = -(gradp(i,j,k,1)+l_gp0[1])*rhoinv + l_gravity[1];
            vel_f(i,j,k,2) = -(gradp(i,j,k,2)+l_gp0[2])*rhoinv + l_gravity[2];
        });
    }
}
