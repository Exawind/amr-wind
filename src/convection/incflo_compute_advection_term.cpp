#include <incflo.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB_utils.H>
#endif

// Need this for convective difference routine
#include "AMReX_MultiFabUtil_3D_C.H"

using namespace amrex;

void
incflo::compute_convective_term (Box const& bx, int lev, MFIter const& mfi,
                                 Array4<Real> const& dvdt, // velocity
                                 Array4<Real> const& drdt, // density
                                 Array4<Real> const& dtdt, // tracer
                                 Array4<Real const> const& vel,
                                 Array4<Real const> const& rho,
                                 Array4<Real const> const& tra,
                                 Array4<Real const> const& umac,
                                 Array4<Real const> const& vmac,
                                 Array4<Real const> const& wmac)
{
    bool regular = true;
#ifdef AMREX_USE_EB
    auto const& fact = EBFactory(lev);
    EBCellFlagFab const& flagfab = fact.getMultiEBCellFlagFab()[mfi];
    Array4<EBCellFlag const> const& flag = flagfab.const_array();
    if (flagfab.getType(bx) == FabType::covered)
    {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            dvdt(i,j,k,0) = 0.0;
            dvdt(i,j,k,1) = 0.0;
            dvdt(i,j,k,2) = 0.0;
            drdt(i,j,k) = 0.0;
        });
        if (advect_tracer) {
            amrex::ParallelFor(bx, ntrac, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                dtdt(i,j,k,n) = 0.0;
            });
        }
        return;
    }

    regular = (flagfab.getType(amrex::grow(bx,1)) == FabType::regular);

    Array4<Real const> fcx, fcy, fcz, vfrac, apx, apy, apz;
    if (!regular) {
        fcx = fact.getFaceCent()[0]->const_array(mfi);
        fcy = fact.getFaceCent()[1]->const_array(mfi);
        fcz = fact.getFaceCent()[2]->const_array(mfi);
        vfrac = fact.getVolFrac().const_array(mfi);
        apx = fact.getAreaFrac()[0]->const_array(mfi);
        apy = fact.getAreaFrac()[1]->const_array(mfi);
        apz = fact.getAreaFrac()[2]->const_array(mfi);
    }
#endif

    int nmaxcomp = AMREX_SPACEDIM;
    if (advect_tracer) nmaxcomp = std::max(nmaxcomp,ntrac);
    Box tmpbox = amrex::surroundingNodes(bx);
    int tmpcomp = nmaxcomp*AMREX_SPACEDIM;
    Box rhotrac_box = amrex::grow(bx,2);
#ifdef AMREX_USE_EB
    Box gbx = bx;
    if (!regular) {
        gbx.grow(2);
        tmpbox.grow(3);
        tmpcomp += nmaxcomp;
        rhotrac_box.grow(2);
    }
#endif

    FArrayBox rhotracfab;
    Elixir eli_rt;
    Array4<Real> rhotrac;
    if (advect_tracer) {
        rhotracfab.resize(rhotrac_box, ntrac);
        eli_rt = rhotracfab.elixir();
        rhotrac = rhotracfab.array();
        amrex::ParallelFor(rhotrac_box, ntrac,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            rhotrac(i,j,k,n) = rho(i,j,k) * tra(i,j,k,n);
        });
    }

    FArrayBox tmpfab(tmpbox, tmpcomp);
    Elixir eli = tmpfab.elixir();
    Array4<Real> fx = tmpfab.array(0);
    Array4<Real> fy = tmpfab.array(nmaxcomp);
    Array4<Real> fz = tmpfab.array(nmaxcomp*2);

#ifdef AMREX_USE_EB
    if (!regular)
    {
        Array4<Real> scratch = tmpfab.array(0);
        Array4<Real> qface = tmpfab.array(nmaxcomp*3);
        Array4<Real> dUdt_tmp = tmpfab.array(nmaxcomp*3);

        // velocity
        compute_convective_fluxes_eb(lev, gbx, AMREX_SPACEDIM, fx, fy, fz, vel, umac, vmac, wmac,
                                     get_velocity_bcrec().data(),
                                     get_velocity_bcrec_device_ptr(),
                                     flag, fcx, fcy, fcz, qface);
        compute_convective_rate_eb(lev, gbx, AMREX_SPACEDIM, dUdt_tmp, fx, fy, fz,
                                   flag, vfrac, apx, apy, apz);
        redistribute_eb(lev, bx, AMREX_SPACEDIM, dvdt, dUdt_tmp, scratch, flag, vfrac);

        // density
        if (!constant_density) {
            compute_convective_fluxes_eb(lev, gbx, 1, fx, fy, fz, rho, umac, vmac, wmac,
                                         get_density_bcrec().data(),
                                         get_density_bcrec_device_ptr(),
                                         flag, fcx, fcy, fcz, qface);
            compute_convective_rate_eb(lev, gbx, 1, dUdt_tmp, fx, fy, fz,
                                       flag, vfrac, apx, apy, apz);
            redistribute_eb(lev, bx, 1, drdt, dUdt_tmp, scratch, flag, vfrac);
        }

        if (advect_tracer) {
            compute_convective_fluxes_eb(lev, gbx, ntrac, fx, fy, fz, rhotrac, umac, vmac, wmac,
                                         get_tracer_bcrec().data(),
                                         get_tracer_bcrec_device_ptr(),
                                         flag, fcx, fcy, fcz, qface);
            compute_convective_rate_eb(lev, gbx, ntrac, dUdt_tmp, fx, fy, fz,
                                       flag, vfrac, apx, apy, apz);
            redistribute_eb(lev, bx, ntrac, dtdt, dUdt_tmp, scratch, flag, vfrac);
        }
    }
    else
#endif
    {
        // velocity
        compute_convective_fluxes(lev, bx, AMREX_SPACEDIM, fx, fy, fz, vel, umac, vmac, wmac,
                                  get_velocity_bcrec().data(),
                                  get_velocity_bcrec_device_ptr());
        compute_convective_rate(lev, bx, AMREX_SPACEDIM, dvdt, fx, fy, fz);

        // density
        if (!constant_density) {
            compute_convective_fluxes(lev, bx, 1, fx, fy, fz, rho, umac, vmac, wmac,
                                      get_density_bcrec().data(),
                                      get_density_bcrec_device_ptr());
            compute_convective_rate(lev, bx, 1, drdt, fx, fy, fz);
        }

        // tracer
        if (advect_tracer) {
            compute_convective_fluxes(lev, bx, ntrac, fx, fy, fz, rhotrac, umac, vmac, wmac,
                                      get_tracer_bcrec().data(),
                                      get_tracer_bcrec_device_ptr());
            compute_convective_rate(lev, bx, ntrac, dtdt, fx, fy, fz);
        }
    }
   
}

void
incflo::compute_convective_term (Vector<MultiFab*> const& conv_u,
                                 Vector<MultiFab*> const& conv_r,
                                 Vector<MultiFab*> const& conv_t,
                                 Vector<MultiFab const*> const& vel,
                                 Vector<MultiFab const*> const& density,
                                 Vector<MultiFab const*> const& tracer,
                                 Vector<MultiFab const*> const& vel_forces,
                                 Vector<MultiFab const*> const& tra_forces,
                                 Real time)
{
    Vector<MultiFab> u_mac(finest_level+1), v_mac(finest_level+1), w_mac(finest_level+1);
    int ngmac = 0;  // This will change for Godunov.
#ifdef AMREX_USE_EB
    if (!EBFactory(0).isAllRegular()) ngmac = 3;
#endif

    for (int lev = 0; lev <= finest_level; ++lev) {
        const BoxArray& ba = density[lev]->boxArray();
        const DistributionMapping& dm = density[lev]->DistributionMap();
        u_mac[lev].define(amrex::convert(ba,IntVect::TheDimensionVector(0)), dm,
                          1, ngmac, MFInfo(), Factory(lev));
        v_mac[lev].define(amrex::convert(ba,IntVect::TheDimensionVector(1)), dm,
                          1, ngmac, MFInfo(), Factory(lev));
        w_mac[lev].define(amrex::convert(ba,IntVect::TheDimensionVector(2)), dm,
                          1, ngmac, MFInfo(), Factory(lev));
        if (ngmac > 0) {
            u_mac[lev].setVal(0.0);
            v_mac[lev].setVal(0.0);
            w_mac[lev].setVal(0.0);
        }
        // Predict normal velocity to faces -- note that the {u_mac, v_mac, w_mac}
        //    returned from this call are on face CENTROIDS
        if (use_godunov) {
            predict_godunov(lev, time, u_mac[lev], v_mac[lev], w_mac[lev], *vel[lev], *density[lev],
                            *vel_forces[lev]);
        } else {
            predict_vels_on_faces(lev, time, u_mac[lev], v_mac[lev], w_mac[lev], *vel[lev]);
        }
    }

    if (use_godunov) {
        VisMF::Write(u_mac[0], "u_mac");
        VisMF::Write(v_mac[0], "v_mac");
        VisMF::Write(w_mac[0], "w_mac");
        amrex::Abort("after predict mac");
    }

    apply_MAC_projection(u_mac, v_mac, w_mac, density, time);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        if (ngmac > 0) {
            u_mac[lev].FillBoundary(geom[lev].periodicity());
            v_mac[lev].FillBoundary(geom[lev].periodicity());
            w_mac[lev].FillBoundary(geom[lev].periodicity());
        }

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling(IntVect(1024,16,16)).SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*density[lev],mfi_info); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            compute_convective_term(bx, lev, mfi,
                                    conv_u[lev]->array(mfi),
                                    conv_r[lev]->array(mfi),
                                    (ntrac>0) ? conv_t[lev]->array(mfi) : Array4<Real>{},
                                    vel[lev]->const_array(mfi),
                                    density[lev]->array(mfi),
                                    (ntrac>0) ? tracer[lev]->array(mfi) : Array4<Real const>{},
                                    u_mac[lev].const_array(mfi),
                                    v_mac[lev].const_array(mfi),
                                    w_mac[lev].const_array(mfi));
        }
    }

}

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

        GpuArray<int,3> iconserv_vel{0,0,0};
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

            GpuArray<int,3> iconserv_density{1,1,1};

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

            GpuArray<int,3> iconserv_trac{1,1,1};
            AMREX_ALWAYS_ASSERT(ntrac <= 3);

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


