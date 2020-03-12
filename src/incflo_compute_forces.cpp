#include <incflo.H>

using namespace amrex;

void incflo::compute_tra_forces (Vector<MultiFab*> const& tra_forces)
{
    // For now we don't have any external forces on the scalars
    if (m_advect_tracer) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            tra_forces[lev]->setVal(0.0);
        }
    }
}

void incflo::compute_vel_forces (Vector<MultiFab*> const& vel_forces,
                                 Vector<MultiFab const*> const& velocity,
                                 Vector<MultiFab const*> const& density,
                                 Vector<MultiFab const*> const& tracer)
{
    for (int lev = 0; lev <= finest_level; ++lev)  
       compute_vel_forces_on_level (lev, *vel_forces[lev], *velocity[lev], *density[lev], *tracer[lev]);
}

void incflo::compute_vel_forces_on_level (int lev,
                                                MultiFab& vel_forces,
                                          const MultiFab& velocity,
                                          const MultiFab& density,
                                          const MultiFab& tracer)
{
    GpuArray<Real,3> l_gravity{m_gravity[0],m_gravity[1],m_gravity[2]};
    GpuArray<Real,3> l_gp0{m_gp0[0], m_gp0[1], m_gp0[2]};
        
    GpuArray<Real,3> east = {m_east[0],m_east[1],m_east[2]};
    GpuArray<Real,3> north = {m_north[0],m_north[1],m_north[2]};
    GpuArray<Real,3> up = {m_up[0],m_up[1],m_up[2]};

    const Real sinphi = m_sinphi;
    const Real cosphi = m_cosphi;
    const Real corfac = m_corfac;

    const Real u = m_ic_u;
    const Real v = m_ic_v;

    const Real umean = m_vx_mean_forcing;//fixme get rid of this global storage
    const Real vmean = m_vy_mean_forcing;//fixme get rid of this global storage

    const Real T0 = m_temperature_values[0];
    const Real thermalExpansionCoeff = m_thermalExpansionCoeff;
    AMREX_ALWAYS_ASSERT(thermalExpansionCoeff>0.0);
        
    const Real dt = m_dt;
    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (int lev = 0; lev <= finest_level; ++lev) {
        for (MFIter mfi(vel_forces,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real> const& vel_f = vel_forces.array(mfi);
            Array4<Real const> const& rho = density.const_array(mfi);
            Array4<Real const> const& gradp = m_leveldata[lev]->gp.const_array(mfi);

            if (m_use_boussinesq) {
                // This uses a Boussinesq approximation where the buoyancy depends on
                //      first tracer rather than density
                Array4<Real const> const& tra = tracer.const_array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real rhoinv = 1.0/rho(i,j,k);
                    Real ft = thermalExpansionCoeff*(T0-tra(i,j,k));
                    vel_f(i,j,k,0) = -gradp(i,j,k,0)*rhoinv + l_gravity[0] * ft;
                    vel_f(i,j,k,1) = -gradp(i,j,k,1)*rhoinv + l_gravity[1] * ft;
                    vel_f(i,j,k,2) = -gradp(i,j,k,2)*rhoinv + l_gravity[2] * ft;
                });
            } else {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real rhoinv = 1.0/rho(i,j,k);
                    vel_f(i,j,k,0) = -(gradp(i,j,k,0)+l_gp0[0])*rhoinv + l_gravity[0];
                    vel_f(i,j,k,1) = -(gradp(i,j,k,1)+l_gp0[1])*rhoinv + l_gravity[1];
                    vel_f(i,j,k,2) = -(gradp(i,j,k,2)+l_gp0[2])*rhoinv + l_gravity[2];
                });
            }
            
            if(m_use_coriolis){
                
                Array4<Real const> const& vel = velocity.const_array(mfi);
                
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    
                    const Real ue = east[0]*vel(i,j,k,0) + east[1]*vel(i,j,k,1) + east[2]*vel(i,j,k,2);
                    const Real un = north[0]*vel(i,j,k,0) + north[1]*vel(i,j,k,1) + north[2]*vel(i,j,k,2);
                    const Real uu = up[0]*vel(i,j,k,0) + up[1]*vel(i,j,k,1) + up[2]*vel(i,j,k,2);

                    const Real ae = +corfac*(un*sinphi - uu*cosphi);
                    const Real an = -corfac*ue*sinphi;
                    const Real au = +corfac*ue*cosphi;

                    const Real ax = ae*east[0] + an*north[0] + au*up[0];
                    const Real ay = ae*east[1] + an*north[1] + au*up[1];
                    const Real az = ae*east[2] + an*north[2] + au*up[2];//fixme do we turn the z component off?
                    
                    vel_f(i,j,k,0) += ax;
                    vel_f(i,j,k,1) += ay;
                    vel_f(i,j,k,2) += az;
                    
                });
            }
            
            if(m_use_abl_forcing){
                const Real dudt = (u-umean)/dt;// fixme make sure this is the correct dt and not dt/2
                const Real dvdt = (v-vmean)/dt;
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    vel_f(i,j,k,0) += dudt;
                    vel_f(i,j,k,1) += dvdt;
                    // vel_f(i,j,k,2) -= 0.0; // assuming periodic in x,y-dir so do not drive flow in z-dir
                });
            }
            
        }
    }
}
