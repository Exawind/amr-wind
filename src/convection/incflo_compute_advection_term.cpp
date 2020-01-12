#include <incflo.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB_utils.H>
#endif

// Need this for convective difference routine
#include "AMReX_MultiFabUtil_3D_C.H"

using namespace amrex;

//
// Compute the three components of the advection term
//
void
incflo::incflo_compute_advection_term( Vector< std::unique_ptr<MultiFab> >& conv_u_in, 
                                       Vector< std::unique_ptr<MultiFab> >& conv_r_in,
                                       Vector< std::unique_ptr<MultiFab> >& conv_t_in,
                                       Vector< std::unique_ptr<MultiFab> >& vel_forces_in,
                                       Vector< std::unique_ptr<MultiFab> >& scal_forces_in,
                                       Vector< std::unique_ptr<MultiFab> >& vel_in,
                                       Vector< std::unique_ptr<MultiFab> >& density_in,
                                       Vector< std::unique_ptr<MultiFab> >& tracer_in,
                                       Real time)
{
    BL_PROFILE("incflo::incflo_compute_advection_term");

    // Temporaries to store fluxes 
    Vector< std::unique_ptr<MultiFab> > fx;
    Vector< std::unique_ptr<MultiFab> > fy;
    Vector< std::unique_ptr<MultiFab> > fz;

    fx.resize(finest_level+1);
    fy.resize(finest_level+1);
    fz.resize(finest_level+1);

    int num_comp;
    int flux_ngrow = 2;

    for (int lev=0; lev <= finest_level; ++lev)
    {
#ifdef AMREX_USE_EB
    const FabFactory<FArrayBox>& factory =  *ebfactory[lev];
#else
    const FabFactory<FArrayBox>& factory = FArrayBoxFactory();
#endif
        if (advect_tracer)
        {
            // Convert tracer to (rho * tracer) so we can use conservative update
            for (int i = 0; i < ntrac; i++) 
               MultiFab::Multiply(*tracer_in[lev],*density_in[lev],0,i,1,tracer_in[lev]->nGrow());
        }

        // **************************************************
        // Compute div (u u) -- the update for velocity
        // **************************************************
        num_comp = 3; 
 
        fx[lev].reset(new MultiFab(m_u_mac[lev]->boxArray(),dmap[lev],num_comp,flux_ngrow,MFInfo(),factory));
        fy[lev].reset(new MultiFab(m_v_mac[lev]->boxArray(),dmap[lev],num_comp,flux_ngrow,MFInfo(),factory));
        fz[lev].reset(new MultiFab(m_w_mac[lev]->boxArray(),dmap[lev],num_comp,flux_ngrow,MFInfo(),factory));

        int iconserv_vel[3] = {0,0,0};
        bool return_state_not_flux;

        
        if (use_godunov) //  This will define the update as (u dot grad u)
            return_state_not_flux = true;
        else //This will define the update as (div(uu))
            return_state_not_flux = false;

        incflo_compute_fluxes(lev, fx, fy, fz, vel_in, 0, vel_forces_in, 0, num_comp,
                              xslopes_u, yslopes_u, zslopes_u, 0,
                              m_u_mac, m_v_mac, m_w_mac, iconserv_vel, return_state_not_flux);

        if (return_state_not_flux) 
            incflo_compute_convective_update  (lev, conv_u_in, fx, fy, fz, num_comp);
        else
            incflo_compute_conservative_update(lev, conv_u_in, fx, fy, fz, num_comp);


        // **************************************************
        // Compute div (rho u) -- the update for density
        // **************************************************
        conv_r_in[lev]->setVal(0.,0,conv_r_in[lev]->nComp(),conv_r_in[lev]->nGrow());

        if (!constant_density)
        {
            num_comp = 1; 

            incflo_compute_slopes(lev, time, *density_in[lev], xslopes_r, yslopes_r, zslopes_r, 0, num_comp);
 
            fx[lev].reset(new MultiFab(m_u_mac[lev]->boxArray(),dmap[lev],num_comp,flux_ngrow,MFInfo(),factory));
            fy[lev].reset(new MultiFab(m_v_mac[lev]->boxArray(),dmap[lev],num_comp,flux_ngrow,MFInfo(),factory));
            fz[lev].reset(new MultiFab(m_w_mac[lev]->boxArray(),dmap[lev],num_comp,flux_ngrow,MFInfo(),factory));

            int iconserv_density[1] = {1};

            // Note that the "ntrac" component of scal_forces holds zeroes
            bool return_state_not_flux = false;
            incflo_compute_fluxes(lev, fx, fy, fz, density_in, 0, scal_forces_in, ntrac, num_comp,
                                  xslopes_r, yslopes_r, zslopes_r, 0,
                                  m_u_mac, m_v_mac, m_w_mac, iconserv_density, return_state_not_flux);

            incflo_compute_conservative_update(lev, conv_r_in, fx, fy, fz, num_comp);

        }

        // **********************************************************
        // Compute div (rho trac u) -- the update for (rho*trac)
        // **********************************************************
        conv_t_in[lev]->setVal(0.,0,conv_t_in[lev]->nComp(),conv_t_in[lev]->nGrow());

        if (advect_tracer)
        {
            num_comp = ntrac; 

            incflo_compute_slopes(lev, time, *tracer_in[lev], xslopes_t, yslopes_t, zslopes_t, 0, num_comp);
 
            fx[lev].reset(new MultiFab(m_u_mac[lev]->boxArray(),dmap[lev],num_comp,flux_ngrow,MFInfo(),factory));
            fy[lev].reset(new MultiFab(m_v_mac[lev]->boxArray(),dmap[lev],num_comp,flux_ngrow,MFInfo(),factory));
            fz[lev].reset(new MultiFab(m_w_mac[lev]->boxArray(),dmap[lev],num_comp,flux_ngrow,MFInfo(),factory));

            int iconserv_trac[ntrac];
            for (int i = 0; i < ntrac; i++) iconserv_trac[i] = 1;

            // Here we are assuming we are updating (rho * trac) conservatively, not trac itself convectively
            bool return_state_not_flux = false;
            incflo_compute_fluxes(lev, fx, fy, fz, tracer_in, 0, scal_forces_in, 0, num_comp,
                                  xslopes_t, yslopes_t, zslopes_t, 0,
                                  m_u_mac, m_v_mac, m_w_mac, iconserv_trac, return_state_not_flux);

            incflo_compute_conservative_update(lev, conv_t_in, fx, fy, fz, num_comp);
        }

        // Convert (rho * tracer) back to tracer
        if (advect_tracer)
        {
            for (int i = 0; i < ntrac; i++) 
               MultiFab::Divide(*tracer_in[lev],*density_in[lev],0,i,1,tracer_in[lev]->nGrow());
        }

        // Return the negative
        conv_u_in[lev] -> mult(-1.0);
        conv_r_in[lev] -> mult(-1.0);
        conv_t_in[lev] -> mult(-1.0);
    } // lev
}

void 
incflo::incflo_compute_conservative_update(const int lev, 
                                           Vector< std::unique_ptr<MultiFab> >& conv_in,
                                           Vector< std::unique_ptr<MultiFab> >& fx,
                                           Vector< std::unique_ptr<MultiFab> >& fy,
                                           Vector< std::unique_ptr<MultiFab> >& fz,
                                           int num_comp)
{
    // Note conv_tmp needs two ghost cells for the redistribution step.
#ifdef AMREX_USE_EB
    MultiFab conv_tmp(grids[lev], dmap[lev], num_comp, 2, MFInfo(), *ebfactory[lev]);
#else
    MultiFab conv_tmp(grids[lev], dmap[lev], num_comp, 2, MFInfo());
#endif
    conv_tmp.setVal(0.);

    Array<MultiFab*,AMREX_SPACEDIM> fluxes;

    fluxes[0] = fx[lev].get();
    fluxes[1] = fy[lev].get();
    fluxes[2] = fz[lev].get();

#ifdef AMREX_USE_EB
    bool already_on_centroids = true;
    EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes), geom[lev], already_on_centroids);
    amrex::single_level_redistribute( lev, conv_tmp, *conv_in[lev], 0, num_comp, geom);
#else
    computeDivergence(*conv_in[lev], GetArrOfConstPtrs(fluxes), geom[lev]);
#endif
}


void 
incflo::incflo_compute_convective_update(const int lev, 
                                         Vector< std::unique_ptr<MultiFab> >& conv,
                                         Vector< std::unique_ptr<MultiFab> >& fx,
                                         Vector< std::unique_ptr<MultiFab> >& fy,
                                         Vector< std::unique_ptr<MultiFab> >& fz,
                                         int num_comp)
{
    const GpuArray<Real,AMREX_SPACEDIM> dxinv = geom[lev].InvCellSizeArray();

    for (MFIter mfi(*conv[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       // Tilebox
       Box bx = mfi.tilebox ();

       Array4<Real const> const& umac_arr = m_u_mac[lev]->const_array(mfi);
       Array4<Real const> const& vmac_arr = m_v_mac[lev]->const_array(mfi);
       Array4<Real const> const& wmac_arr = m_w_mac[lev]->const_array(mfi);

       Array4<Real const> const& fx_arr = fx[lev]->const_array(mfi);
       Array4<Real const> const& fy_arr = fy[lev]->const_array(mfi);
       Array4<Real const> const& fz_arr = fz[lev]->const_array(mfi);

       Array4<Real> const& diff_arr = conv[lev]->array(mfi);

       amrex_compute_convective_difference(bx, diff_arr, umac_arr, vmac_arr, wmac_arr, fx_arr, fy_arr, fz_arr, dxinv);
    }
}


