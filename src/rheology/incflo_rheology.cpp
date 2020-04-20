#include <incflo.H>
#include <derive_K.H>

using namespace amrex;

namespace {

struct NonNewtonianEddyViscosity
{
    FluidModel fluid_model;
    amrex::Real mu, smag_const;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (amrex::Real sr, amrex::Real den, amrex::Real ds) const noexcept {
        switch (fluid_model)
        {
        case FluidModel::SmagorinskyLillySGS:
        {
            // fixme this is not good place for Smagorinsky move to an ABL specific location
            return den*smag_const*smag_const*ds*ds*sr;
        }
        default:
        {
            return mu;
        }
        };
    }
};

}

void incflo::compute_viscosity (Vector<MultiFab*> const& vel_eta,
                                Vector<MultiFab*> const& tra_eta,
                                Vector<MultiFab const*> const& rho,
                                Vector<MultiFab const*> const& vel,
                                Vector<MultiFab const*> const& /* tra */,
                                Real /* time */, int nghost)
{
    if (m_fluid_model == FluidModel::Newtonian)
    {
        for (auto mf : vel_eta) {
            mf->setVal(m_mu, 0, 1, nghost);
        }
    }
    else
    {
        NonNewtonianEddyViscosity non_newtonian_eddy_viscosity;
        non_newtonian_eddy_viscosity.fluid_model = m_fluid_model;
        non_newtonian_eddy_viscosity.mu = m_mu;
        non_newtonian_eddy_viscosity.smag_const = m_Smagorinsky_Lilly_SGS_constant;

        for (int lev = 0; lev <= finest_level; ++lev) {

            const Real dx = geom[lev].CellSize()[0];
            const Real dy = geom[lev].CellSize()[1];
            const Real dz = geom[lev].CellSize()[2];
            const Real ds = pow(dx*dy*dz,1.0/3.0);

            Real idx = 1.0 / dx;
            Real idy = 1.0 / dy;
            Real idz = 1.0 / dz;
            
            const Geometry& gm = Geom(lev);
            const Box& domain = gm.Domain();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*vel_eta[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.growntilebox(nghost);
                Array4<Real> const& eta_arr = vel_eta[lev]->array(mfi);
                Array4<Real const> const& vel_arr = vel[lev]->const_array(mfi);
                Array4<Real const> const& rho_arr = rho[lev]->const_array(mfi);

                
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real sr = incflo_strainrate<StencilInterior>(i,j,k,idx,idy,idz,vel_arr);
                    Real den = rho_arr(i,j,k);
                    eta_arr(i,j,k) = non_newtonian_eddy_viscosity(sr,den,ds);

                });
                
                
                /* one-sided stencils for non-periodic bc's */                
                Box const& bxi = mfi.tilebox();

                int idim = 0;
                if (!gm.isPeriodic(idim)) {
                    if (bxi.smallEnd(idim) == domain.smallEnd(idim)) {
                        
                        IntVect low(bxi.smallEnd());
                        IntVect hi(bxi.bigEnd());
                        int sm = low[idim];
                        low.setVal(idim,sm);
                        hi.setVal(idim,sm);
                      
                        Box bxlo = Box(low,hi).grow({0,1,1});

                        amrex::ParallelFor(bxlo,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            Real sr = incflo_strainrate<StencilILO>(i,j,k,idx,idy,idz,vel_arr);
                            Real den = rho_arr(i,j,k);
                            eta_arr(i,j,k) = non_newtonian_eddy_viscosity(sr,den,ds);
                        });
                    }
                    
                    if (bxi.bigEnd(idim) == domain.bigEnd(idim)) {
                        
                        IntVect low(bxi.smallEnd());
                        IntVect hi(bxi.bigEnd());
                        int sm = hi[idim];
                        low.setVal(idim,sm);
                        hi.setVal(idim,sm);

                        Box bxhi = Box(low,hi).grow({0,1,1});

                        amrex::ParallelFor(bxhi,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                           Real sr = incflo_strainrate<StencilIHI>(i,j,k,idx,idy,idz,vel_arr);
                           Real den = rho_arr(i,j,k);
                           eta_arr(i,j,k) = non_newtonian_eddy_viscosity(sr,den,ds);
                        });
                    }
                }
                
                idim = 1;
                if (!gm.isPeriodic(idim)) {
                    if (bxi.smallEnd(idim) == domain.smallEnd(idim)) {
                        
                        IntVect low(bxi.smallEnd());
                        IntVect hi(bxi.bigEnd());
                        int sm = low[idim];
                        low.setVal(idim,sm);
                        hi.setVal(idim,sm);
                      
                        Box bxlo = Box(low,hi).grow({1,0,1});

                        amrex::ParallelFor(bxlo,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            Real sr = incflo_strainrate<StencilJLO>(i,j,k,idx,idy,idz,vel_arr);
                            Real den = rho_arr(i,j,k);
                            eta_arr(i,j,k) = non_newtonian_eddy_viscosity(sr,den,ds);
                        });
                    }
                    
                    if (bxi.bigEnd(idim) == domain.bigEnd(idim)) {
                        
                        IntVect low(bxi.smallEnd());
                        IntVect hi(bxi.bigEnd());
                        int sm = hi[idim];
                        low.setVal(idim,sm);
                        hi.setVal(idim,sm);

                        Box bxhi = Box(low,hi).grow({1,0,1});

                        amrex::ParallelFor(bxhi,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                           Real sr = incflo_strainrate<StencilJHI>(i,j,k,idx,idy,idz,vel_arr);
                           Real den = rho_arr(i,j,k);
                           eta_arr(i,j,k) = non_newtonian_eddy_viscosity(sr,den,ds);
                        });
                    }
                }
                
                idim = 2;
                if (!gm.isPeriodic(idim)) {
                    if (bxi.smallEnd(idim) == domain.smallEnd(idim)) {
                        
                        IntVect low(bxi.smallEnd());
                        IntVect hi(bxi.bigEnd());
                        int sm = low[idim];
                        low.setVal(idim,sm);
                        hi.setVal(idim,sm);
                      
                        Box bxlo = Box(low,hi).grow({1,1,0});

                        amrex::ParallelFor(bxlo,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            Real sr = incflo_strainrate<StencilKLO>(i,j,k,idx,idy,idz,vel_arr);
                            Real den = rho_arr(i,j,k);
                            eta_arr(i,j,k) = non_newtonian_eddy_viscosity(sr,den,ds);
                        });
                    }
                    
                    if (bxi.bigEnd(idim) == domain.bigEnd(idim)) {
                        
                        IntVect low(bxi.smallEnd());
                        IntVect hi(bxi.bigEnd());
                        int sm = hi[idim];
                        low.setVal(idim,sm);
                        hi.setVal(idim,sm);

                        Box bxhi = Box(low,hi).grow({1,1,0});

                        amrex::ParallelFor(bxhi,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                           Real sr = incflo_strainrate<StencilKHI>(i,j,k,idx,idy,idz,vel_arr);
                           Real den = rho_arr(i,j,k);
                           eta_arr(i,j,k) = non_newtonian_eddy_viscosity(sr,den,ds);
                        });
                    }
                }
            }
        }
    }

    switch(m_fluid_model){
        case FluidModel::SmagorinskyLillySGS:
        {

            for (auto mf : tra_eta) {
                for (int n = 0; n < m_ntrac; ++n) {
                    mf->setVal(0.0, n, 1, nghost);
                }
            }

            if(m_ntrac && m_advect_tracer){

                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_ntrac == 1,"SmagorinskyLillySGS only implemented for 1 tracer");


                const Real Pr = 0.70; // fixme make an input
                const Real Pr_t = 0.3333; // fixme make an input

                for (int lev = 0; lev <= finest_level; ++lev) {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                    for (MFIter mfi(*vel_eta[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                    {
                        Box const& bx = mfi.growntilebox(nghost);
                        Array4<Real> const& vel_eta_arr = vel_eta[lev]->array(mfi);
                        Array4<Real> const& tra_eta_arr = tra_eta[lev]->array(mfi);

                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {

                            const Real mu_t = vel_eta_arr(i,j,k);
                            tra_eta_arr(i,j,k) = m_mu/Pr + mu_t/Pr_t;
                            vel_eta_arr(i,j,k) = m_mu + mu_t;

                        });
                    }
                }
            }


            break;
        }
        default:
        {
            for (auto mf : tra_eta) {
                for (int n = 0; n < m_ntrac; ++n) {
                    mf->setVal(m_mu_s[n], n, 1, nghost);
                }
            }
        }
    };
   

}
