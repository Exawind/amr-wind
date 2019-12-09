#include "incflo.H"
#include <random>
#include <boundary_conditions_F.H>

void
incflo::init_abl(MultiFab& density_mfab, MultiFab& vel_mfab, MultiFab& tracer_mfab, Real dx, Real dy, Real dz)
{
    const Real pi = std::acos(-1.0);
    const Real aval = Uperiods*2.0*pi/(geom[0].ProbHi(1) - geom[0].ProbLo(1));
    const Real bval = Vperiods*2.0*pi/(geom[0].ProbHi(0) - geom[0].ProbLo(0));
    const Real ufac = deltaU*std::exp(0.5)/zRefHeight;
    const Real vfac = deltaV*std::exp(0.5)/zRefHeight;
    
    AMREX_ASSERT(ntrac>0);
    
    for (MFIter mfi(vel_mfab,false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        
        const auto& den_arr = density_mfab.array(mfi);
        const auto& vel_arr = vel_mfab.array(mfi);
        const auto& trac_arr = tracer_mfab.array(mfi);

        for(int k(bx.loVect()[2]); k <= bx.hiVect()[2]; k++){
            for(int j(bx.loVect()[1]); j <= bx.hiVect()[1]; j++){
                for(int i(bx.loVect()[0]); i <= bx.hiVect()[0]; i++){
    
                    // constant density
                    den_arr(i,j,k) = ro_0;

                    // physical location of cell center
                    const Real x = geom[0].ProbLo(0) + (i+0.5)*dx;
                    const Real y = geom[0].ProbLo(1) + (j+0.5)*dy;
                    const Real z = geom[0].ProbLo(2) + (k+0.5)*dz;

                    const Real zl = z/zRefHeight;
                    const Real damp = std::exp(-0.5*zl*zl);

                    vel_arr(i,j,k,0) = ic_u;
                    vel_arr(i,j,k,1) = ic_v;
                    vel_arr(i,j,k,2) = ic_w;

                    // velocity field plus perturbations
                    if(z < fabs(ic_u)) vel_arr(i,j,k,0) = copysign(z,ic_u);
                    if(z < fabs(ic_v)) vel_arr(i,j,k,1) = copysign(z,ic_v);

                    vel_arr(i,j,k,0) += ufac*damp*z*std::cos(aval*y);
                    vel_arr(i,j,k,1) += vfac*damp*z*std::sin(bval*x);

                    // potential temperature and boussinesq stuff
                    Real theta = temperature_values[0];
                    for(int t = 0; t < ntemperature-1; ++t){
                        const Real slope = (temperature_values[t+1]-temperature_values[t])/(temperature_heights[t+1]-temperature_heights[t]);
                        if(z > temperature_heights[t] && z <= temperature_heights[t+1]){
                          theta = temperature_values[t] + (z-temperature_heights[t])*slope;
                        }
                    }

                    // add perturbations to potential temperature
        //            if(z < cutoff_height){
        //                // Random number generator for adding temperature perturbations
        //                const Real thetaGaussMean_ = 0.0;
        //                const Real thetaGaussVar_ = 1.0;
        //                std::random_device rd{};
        //                std::mt19937 gen{rd()};
        //                std::normal_distribution<Real> gaussNormal(thetaGaussMean_, thetaGaussVar_);
        //                theta += theta_amplitude*gaussNormal(gen);
        //            }

                    // initialize all tracers to be zero
                    for(int n=0;n<ntrac;++n){
                        trac_arr(i,j,k,n)  = 0.0;
                    }

                    // initialize first tracer to potential density
                    if(ntrac>0){
                        trac_arr(i,j,k,0) = theta;
                    }

                    // initialize second tracer to total kinetic energy for 1 eqn SGS model
                    if(ntrac>1){
                        trac_arr(i,j,k,1) = 0.5*(pow(vel_arr(i,j,k,0),2) + pow(vel_arr(i,j,k,1),2) + pow(vel_arr(i,j,k,2),2));
                    }

                } // i loop
            } // j loop
        } // k loop
    } // mfiter loop
}

// has support for filling in ghosts but probably not needed?
void
incflo::add_abl_source_terms(MultiFab& vel_in, MultiFab& tracer_in, int nghost)
{

    AMREX_ASSERT(vel_in.nGrow() >= nghost);
    AMREX_ASSERT(vel_in.nGrow() == tracer_in.nGrow() );
    AMREX_ASSERT(ntrac > 0);
    
    BL_PROFILE("incflo::add_abl_source_terms()");

    // these could all be combined in to one function if we wanted...
    
    if(use_boussinesq)  add_boussinesq(vel_in, tracer_in, nghost);
    if(coriolis_effect) add_coriolis(vel_in,nghost);
    if(abl_forcing)     add_abl_forcing(vel_in, nghost);
    
}


void
incflo::add_boussinesq(MultiFab& vel_in, MultiFab& tracer_in, int nghost)
{

    AMREX_ASSERT(vel_in.nGrow() >= nghost);
    AMREX_ASSERT(vel_in.nGrow() == tracer_in.nGrow() );

    BL_PROFILE("incflo::add_boussinesq()");

    const auto dt_ = dt;
    const auto T0 = temperature_values[0];
    const auto thermalExpansionCoeff_ = thermalExpansionCoeff;
    amrex::GpuArray<Real,3> g;
    g[0] = gravity[0];
    g[1] = gravity[1];
    g[2] = gravity[2];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(vel_in,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
     
        const auto v = vel_in.array(mfi);
        const auto T = tracer_in.array(mfi);
        
        amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real fact = dt_*thermalExpansionCoeff_*(T0 - T(i,j,k,0));

            v(i,j,k,0) += fact*g[0];
            v(i,j,k,1) += fact*g[1];
            v(i,j,k,2) += fact*g[2];
            
         });
    }
}


void
incflo::add_coriolis(MultiFab& vel_in, int nghost)
{

    AMREX_ASSERT(vel_in.nGrow() >= nghost);

    BL_PROFILE("incflo::add_coriolis()");

    // fixme could move somewhere global since this is a constant
    const Real sinphi = std::sin(latitude);
    const Real cosphi = std::cos(latitude);
    const Real dt_ = dt;
    const Real corfac_ = corfac;
    amrex::GpuArray<Real,3> e, n, u;
    e = {1.0,0.0,0.0};
    n = {0.0,1.0,0.0};
    u = {0.0,0.0,1.0};
    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(vel_in,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
     
        auto const v = vel_in.array(mfi);
        
        amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
                            
            const Real ue = e[0]*v(i,j,k,0) + e[1]*v(i,j,k,1) + e[2]*v(i,j,k,2);
            const Real un = n[0]*v(i,j,k,0) + n[1]*v(i,j,k,1) + n[2]*v(i,j,k,2);
            const Real uu = u[0]*v(i,j,k,0) + u[1]*v(i,j,k,1) + u[2]*v(i,j,k,2);

            const Real ae = +corfac_*( un*sinphi - uu*cosphi);
            const Real an = -corfac_*ue*sinphi;
            const Real au = +corfac_*ue*cosphi;

            const Real ax = ae*e[0] + an*n[0] + au*u[0];
            const Real ay = ae*e[1] + an*n[1] + au*u[1];
            const Real az = ae*e[2] + an*n[2] + au*u[2];

            // add in coriolis without density since that is done in the next step
            v(i,j,k,0) -= dt_*ax;
            v(i,j,k,1) -= dt_*ay;
            v(i,j,k,2) -= 0*dt_*az;//fixme there might be known issues with using z component turn off for now
            
        });
    }
}

void
incflo::add_abl_forcing(MultiFab& vel_in, int nghost)
{
    AMREX_ASSERT(vel_in.nGrow() >= nghost);

    BL_PROFILE("incflo::add_abl_forcing()");

    const Real ic_u_ = ic_u;
    const Real ic_v_ = ic_v;
    const Real vx_mean_ = vx_mean;
    const Real vy_mean_ = vy_mean;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(vel_in,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
     
        const auto v = vel_in.array(mfi);
        
        amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // rho is handled on the outside
            // dt*omega*[ rho*(u-<u>)/dt ]  dt/dt cancel out
            const Real omega = 1.0; // under-relaxation
            v(i,j,k,0) += omega*(ic_u_ - vx_mean_);
            v(i,j,k,1) += omega*(ic_v_ - vy_mean_);
            
        });
    }
}



// note this adds to momentum and sets viscosity in tracer
void incflo::add_eddy_viscosity(Vector<std::unique_ptr<MultiFab>>& eta_out,
                                Vector<std::unique_ptr<MultiFab>>& eta_tracer_out,
                                Real time_in)
{
    BL_PROFILE("incflo::add_eddy_viscosity()");

    if(!sgs_model) return;

    if(sgs_model && ntrac == 1) ComputeStrainrate(time_in);
       
    if(ntrac > 1)
    {
        amrex::Abort("ntrac > 1 is for ksgs and not finished yet\n");
    }
    
    const Real zlo = geom[0].ProbLo(2);
    const Real utau_o_nu = utau/(mu/ro_0);
    const Real Prandtl_turb = 0.333333; //fixme make an input
    const Real Ri_critical = 0.3; // 0.2 <= Ric <= 0.4 critical Richardson number fixme should be an input

    for(int lev = 0; lev <= finest_level; lev++)
    {
        Box domain(geom[lev].Domain());
        
        const Real dx = geom[lev].CellSize()[0];
        const Real dy = geom[lev].CellSize()[1];
        const Real dz = geom[lev].CellSize()[2];
        const Real ds = pow(dx*dy*dz,1.0/3.0);

        const Real Cs_ds2 = pow(Smagorinsky_Lilly_SGS_constant*ds,2);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for(MFIter mfi(*eta_out[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox();

            const auto viscosity_arr = eta_out[lev]->array(mfi);
            const auto viscosity_tracer_arr = eta_tracer_out[lev]->array(mfi);
            
            const auto vel_arr = vel[lev]->array(mfi);
            const auto strainrate_arr = strainrate[lev]->array(mfi);
            const auto tracer_arr = tracer[lev]->array(mfi);
            const auto den_arr = density[lev]->array(mfi);
            
            // Smagorinsky–Lilly SGS model
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const Real z = zlo + (k+0.5)*dz;
                const Real yp = z*utau_o_nu;
                const Real damp = pow(1.0-exp(-yp/26.0),2);// this probably does not turn on for ABL

                // make sure ghosts exist before accessing k+1
//              const Real dthetadz = (tracer_arr(i,j,k+1,0)-tracer_arr(i,j,k,0))/dz;
//              const Real dudz = (vel_arr(i,j,k+1,0)-vel_arr(i,j,k,0))/dz;
//              const Real dvdz = (vel_arr(i,j,k+1,1)-vel_arr(i,j,k,1))/dz;
                const Real Ri = 0.0;// = fabs(gravity[2]*dthetadz)/(pow(dudz,2)+pow(dvdz,2)+1.0e-12)/ro_0;
                
                const Real KM = Cs_ds2*pow(1.0-Ri/Ri_critical,0.5)*strainrate_arr(i,j,k);
                const Real KH = KM/Prandtl_turb;
            
                // add eddy viscosity to eta (mu+mu_t) where mu = rho*nu
                viscosity_arr(i,j,k) += KM*damp*den_arr(i,j,k);
                // set eddy viscosity for potential temperature equation
                viscosity_tracer_arr(i,j,k,0) = KH*damp*den_arr(i,j,k);
                
            });
     
        }

        
        // fixme not sure why this is done yet
        eta_out[lev]->FillBoundary(geom[lev].periodicity());
        eta_tracer_out[lev]->FillBoundary(geom[lev].periodicity());
        

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for(MFIter mfi(*density[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            fill_bc0(BL_TO_FORTRAN_ANYD((*eta_out[lev])[mfi]),
                     bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                     bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                     bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                     domain.loVect(), domain.hiVect(),
                     &nghost);
            
            fill_bc0(BL_TO_FORTRAN_ANYD((*eta_tracer_out[lev])[mfi]),
                     bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                     bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                     bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                     domain.loVect(), domain.hiVect(),
                     &nghost);
            
        }

        eta_out[lev]->FillBoundary(geom[lev].periodicity());
        eta_tracer_out[lev]->FillBoundary(geom[lev].periodicity());
    }

}
