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

        AMREX_FOR_3D(bx, i, j, k,
        {
            
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
            if(z < ic_u) vel_arr(i,j,k,0) = z;
            if(z < ic_v) vel_arr(i,j,k,1) = z;

            vel_arr(i,j,k,0) += ufac*damp*z*std::cos(aval*y);
            vel_arr(i,j,k,1) += vfac*damp*z*std::sin(bval*x);
       
            

            // potential temperature and boussinesq stuff
            Real theta = temperature_values[0];
            for(int i = 0; i < ntemperature-1; ++i){
                const Real slope = (temperature_values[i+1]-temperature_values[i])/(temperature_heights[i+1]-temperature_heights[i]);
                if(z > temperature_heights[i] && z <= temperature_heights[i+1]){
                  theta = temperature_values[i] + (z-temperature_heights[i])*slope;
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
            
            
        });
    }
}

// has support for filling in ghosts but probably not needed?
void
incflo::add_abl_source_terms(MultiFab& vel_in, MultiFab& tracer_in, int nghost)
{

    AMREX_ASSERT(vel_in.nGrow() >= nghost);
    AMREX_ASSERT(vel_in.nGrow() == tracer_in.ngrow() );
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
    AMREX_ASSERT(vel_in.nGrow() == tracer_in.ngrow() );

    BL_PROFILE("incflo::add_boussinesq()");

    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(vel_in,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
     
        auto const v = vel_in.array(mfi);
        auto const T = tracer_in.array(mfi);
        
        AMREX_FOR_3D ( bx, i, j, k,
        {
            const Real fact = dt*thermalExpansionCoeff*(temperature_values[0] - T(i,j,k,0));

            v(i,j,k,0) += fact*gravity[0];
            v(i,j,k,1) += fact*gravity[1];
            v(i,j,k,2) += fact*gravity[2];
            
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
	
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(vel_in,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
     
        auto const v = vel_in.array(mfi);
        
        AMREX_FOR_3D ( bx, i, j, k,
        {
                            
            const Real ue = east[0]*v(i,j,k,0)  + east[1]*v(i,j,k,1)  + east[2]*v(i,j,k,2);
            const Real un = north[0]*v(i,j,k,0) + north[1]*v(i,j,k,1) + north[2]*v(i,j,k,2);
            const Real uu = up[0]*v(i,j,k,0)    + up[1]*v(i,j,k,1)    + up[2]*v(i,j,k,2);

            const Real ae = +corfac*( un*sinphi - uu*cosphi);
            const Real an = -corfac*ue*sinphi;
            const Real au = +corfac*ue*cosphi;

            const Real ax = ae*east[0] + an*north[0] + au*up[0];
            const Real ay = ae*east[1] + an*north[1] + au*up[1];
            const Real az = ae*east[2] + an*north[2] + au*up[2];

            // add in coriolis without density since that is done in the next step
            //fixme check the sign
            v(i,j,k,0) -= dt*ax;
            v(i,j,k,1) -= dt*ay;
            v(i,j,k,2) -= 0*dt*az;//fixme there might be known issues with using z component turn off for now
            
        });
    }
}

void
incflo::add_abl_forcing(MultiFab& vel_in, int nghost)
{
    AMREX_ASSERT(vel_in.nGrow() >= nghost);

    BL_PROFILE("incflo::add_abl_forcing()");

  
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(vel_in,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
     
        auto const v = vel_in.array(mfi);
        
        AMREX_FOR_3D ( bx, i, j, k,
        {
            // rho is handled on the outside
            // dt*omega*[ rho*(u-<u>)/dt ]  dt/dt cancel out
            const Real omega = 1.0; // under-relaxation
            v(i,j,k,0) += omega*(ic_u - vx_mean);
            v(i,j,k,1) += omega*(ic_v - vy_mean);
            
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
       
    for(int lev = 0; lev <= finest_level; lev++)
    {
        Box domain(geom[lev].Domain());
        
        const Real dx = geom[lev].CellSize()[0];
        const Real dy = geom[lev].CellSize()[1];
        const Real dz = geom[lev].CellSize()[2];
        const Real ds = pow(dx*dy*dz,1.0/3.0);
        
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for(MFIter mfi(*eta_out[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox();

            const auto& viscosity_arr = eta_out[lev]->array(mfi);
            const auto& viscosity_tracer_arr = eta_tracer_out[lev]->array(mfi);
            
            const auto& vel_arr = vel[lev]->array(mfi);
            const auto& strainrate_arr = strainrate[lev]->array(mfi);
            const auto& tracer_arr = tracer[lev]->array(mfi);
            
            if(ntrac == 1){
                
                // Smagorinsky–Lilly SGS model
                AMREX_FOR_3D ( bx, i, j, k,
                {
                    const Real z = geom[0].ProbLo(2) + (k+0.5)*dz;
                    const Real yp = z*utau/(mu/ro_0);
                    const Real damp = pow(1.0-exp(-0.04*yp),2);

                    const Real Cs = 0.18; // 0.18 <= Cs <= 0.25 fixme should be an input

//                    const Real dthetadz = (tracer_arr(i,j,k+1,0)-tracer_arr(i,j,k,0))/dz;
//                    const Real dudz = (vel_arr(i,j,k+1,0)-vel_arr(i,j,k,0))/dz;
//                    const Real dvdz = (vel_arr(i,j,k+1,1)-vel_arr(i,j,k,1))/dz;
                    const Real Ri = 0.0;//fabs(gravity[2]*dthetadz)/(pow(dudz,2)+pow(dvdz,2)+1.0e-12)/ro_0;
                    const Real Ric = 0.3; // 0.2 <= Ric <= 0.4 fixme should be an input
//                    printf("%f Ri number %f\n",z,Ri);
                    
                    const Real KM = pow(Cs*ds,2)*pow(1.0-Ri/Ric,0.5)*strainrate_arr(i,j,k);
                    const Real Prandtl_turb = 0.333333;
                    const Real KH = KM/Prandtl_turb;
                    
                    viscosity_arr(i,j,k) += KM*damp;
                    viscosity_tracer_arr(i,j,k,0) = KH*damp;
                    
                });
                
            } else if(ntrac > 1)
            {
                amrex::Abort("not finished yet\n");
                // Deardorff TKE model SGS model
                AMREX_FOR_3D ( bx, i, j, k,
                {
                    const Real z = geom[0].ProbLo(2) + (k+0.5)*dz;
                    const Real yp = z*utau/(mu/ro_0);
                    const Real damp = pow(1.0-exp(-0.04*yp),2);
                    
                    const Real l = ds;//fixme need to adjust l based on stratification see Moeng 1984
                    const Real KM = 0.1*l*pow(tracer_arr(i,j,k,1),0.5);
                    const Real KH = (1.0+2.0*l/ds)*KM;
                    
                    viscosity_arr(i,j,k) += damp*KM;
                    viscosity_tracer_arr(i,j,k,0) = damp*KH;
                    viscosity_tracer_arr(i,j,k,1) = 0.0;

                });
                
            }

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
