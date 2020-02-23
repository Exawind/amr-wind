#include <incflo.H>

using namespace amrex;

void incflo::init_advection ()
{
    m_iconserv_velocity.resize(AMREX_SPACEDIM, 0);
    m_iconserv_velocity_d.resize(AMREX_SPACEDIM, 0);

    m_iconserv_density.resize(1, 1);
    m_iconserv_density_d.resize(1, 1);

    m_iconserv_tracer.resize(m_ntrac, 1);
    m_iconserv_tracer_d.resize(m_ntrac, 1);
}

void
incflo::compute_convective_term (Vector<MultiFab*> const& conv_u,
                                 Vector<MultiFab*> const& conv_r,
                                 Vector<MultiFab*> const& conv_t,
                                 Vector<MultiFab const*> const& vel,
                                 Vector<MultiFab const*> const& density,
                                 Vector<MultiFab const*> const& tracer,
                                 Vector<MultiFab*> const& u_mac,
                                 Vector<MultiFab*> const& v_mac,
                                 Vector<MultiFab*> const& w_mac,
                                 Vector<MultiFab const*> const& vel_forces,
                                 Vector<MultiFab const*> const& tra_forces,
                                 Real time)
{
//  Vector<MultiFab> u_mac(finest_level+1), v_mac(finest_level+1), w_mac(finest_level+1);
    int ngmac = nghost_mac();

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Predict normal velocity to faces -- note that the {u_mac, v_mac, w_mac}
        //    returned from this call are on face CENTROIDS
        if (m_use_godunov) {
            predict_godunov(lev, time, *u_mac[lev], *v_mac[lev], *w_mac[lev], *vel[lev], *density[lev],
                            *vel_forces[lev]);
        } else {
            predict_vels_on_faces(lev, time, *u_mac[lev], *v_mac[lev], *w_mac[lev], *vel[lev]);
        }
    }

    apply_MAC_projection(u_mac, v_mac, w_mac, density, time);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        if (ngmac > 0) {
            u_mac[lev]->FillBoundary(geom[lev].periodicity());
            v_mac[lev]->FillBoundary(geom[lev].periodicity());
            w_mac[lev]->FillBoundary(geom[lev].periodicity());
        }

        MFItInfo mfi_info;
        // if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling(IntVect(1024,16,16)).SetDynamic(true);
        if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling(IntVect(1024,1024,1024)).SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*density[lev],mfi_info); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            compute_convective_term(bx, lev, mfi,
                                    conv_u[lev]->array(mfi),
                                    conv_r[lev]->array(mfi),
                                    (m_ntrac>0) ? conv_t[lev]->array(mfi) : Array4<Real>{},
                                    vel[lev]->const_array(mfi),
                                    density[lev]->array(mfi),
                                    (m_ntrac>0) ? tracer[lev]->const_array(mfi) : Array4<Real const>{},
                                    u_mac[lev]->const_array(mfi),
                                    v_mac[lev]->const_array(mfi),
                                    w_mac[lev]->const_array(mfi),
                                    (!vel_forces.empty()) ? vel_forces[lev]->const_array(mfi)
                                                          : Array4<Real const>{},
                                    (!tra_forces.empty()) ? tra_forces[lev]->const_array(mfi)
                                                          : Array4<Real const>{});
        }
    }
}

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
                                 Array4<Real const> const& wmac,
                                 Array4<Real const> const& fvel,
                                 Array4<Real const> const& ftra)
{
#ifdef AMREX_USE_EB
    AMREX_ALWAYS_ASSERT(!m_use_godunov);

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
        if (m_advect_tracer) {
            amrex::ParallelFor(bx, m_ntrac, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                dtdt(i,j,k,n) = 0.0;
            });
        }
        return;
    }

    bool regular = (flagfab.getType(amrex::grow(bx,1)) == FabType::regular);

    Array4<Real const> fcx, fcy, fcz, ccc, vfrac, apx, apy, apz;
    if (!regular) {
        fcx = fact.getFaceCent()[0]->const_array(mfi);
        fcy = fact.getFaceCent()[1]->const_array(mfi);
        fcz = fact.getFaceCent()[2]->const_array(mfi);
        ccc = fact.getCentroid().const_array(mfi);
        vfrac = fact.getVolFrac().const_array(mfi);
        apx = fact.getAreaFrac()[0]->const_array(mfi);
        apy = fact.getAreaFrac()[1]->const_array(mfi);
        apz = fact.getAreaFrac()[2]->const_array(mfi);
    }
#endif

    Box rhotrac_box = amrex::grow(bx,2);
    if (m_use_godunov) rhotrac_box.grow(1);
#ifdef AMREX_USE_EB
    if (!regular) rhotrac_box.grow(2);
#endif

    FArrayBox rhotracfab;
    Elixir eli_rt;
    Array4<Real> rhotrac;
    if (m_advect_tracer) {
        rhotracfab.resize(rhotrac_box, m_ntrac);
        if (!m_use_godunov) {
            eli_rt = rhotracfab.elixir();
        }
        rhotrac = rhotracfab.array();
        amrex::ParallelFor(rhotrac_box, m_ntrac,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            rhotrac(i,j,k,n) = rho(i,j,k) * tra(i,j,k,n);
        });
    }

    int nmaxcomp = AMREX_SPACEDIM;
    if (m_advect_tracer) nmaxcomp = std::max(nmaxcomp,m_ntrac);

    if (m_use_godunov)
    {
        FArrayBox tmpfab(amrex::grow(bx,1), nmaxcomp*14+1);
//        Elixir eli = tmpfab.elixir();

        compute_godunov_advection(lev, bx, AMREX_SPACEDIM,
                                  dvdt, vel,
                                  umac, vmac, wmac, fvel,
                                  get_velocity_bcrec_device_ptr(),
                                  get_velocity_iconserv_device_ptr(),
                                  tmpfab.dataPtr());
        if (!m_constant_density) {
            compute_godunov_advection(lev, bx, 1,
                                      drdt, rho,
                                      umac, vmac, wmac, {},
                                      get_density_bcrec_device_ptr(),
                                      get_density_iconserv_device_ptr(),
                                      tmpfab.dataPtr());
        }
        if (m_advect_tracer) {
            compute_godunov_advection(lev, bx, m_ntrac,
                                      dtdt, tra,
                                      umac, vmac, wmac, ftra,
                                      get_tracer_bcrec_device_ptr(),
                                      get_tracer_iconserv_device_ptr(),
                                      tmpfab.dataPtr());
        }
        Gpu::streamSynchronize();
    }
    else
    {
        Box tmpbox = amrex::surroundingNodes(bx);
        int tmpcomp = nmaxcomp*AMREX_SPACEDIM;
#ifdef AMREX_USE_EB
        Box gbx = bx;
        if (!regular) {
            gbx.grow(2);
            tmpbox.grow(3);
            tmpcomp += nmaxcomp;
        }
#endif

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
            compute_convective_fluxes_eb(lev, gbx, AMREX_SPACEDIM,
                                         fx, fy, fz, vel, umac, vmac, wmac,
                                         get_velocity_bcrec().data(),
                                         get_velocity_bcrec_device_ptr(),
                                         flag, fcx, fcy, fcz, ccc);
            compute_convective_rate_eb(lev, gbx, AMREX_SPACEDIM, dUdt_tmp, fx, fy, fz,
                                       flag, vfrac, apx, apy, apz);
            redistribute_eb(lev, bx, AMREX_SPACEDIM, dvdt, dUdt_tmp, scratch, flag, vfrac);

            // density
            if (!m_constant_density) {
                compute_convective_fluxes_eb(lev, gbx, 1,
                                             fx, fy, fz, rho, umac, vmac, wmac,
                                             get_density_bcrec().data(),
                                             get_density_bcrec_device_ptr(),
                                             flag, fcx, fcy, fcz, ccc);
                compute_convective_rate_eb(lev, gbx, 1, dUdt_tmp, fx, fy, fz,
                                           flag, vfrac, apx, apy, apz);
                redistribute_eb(lev, bx, 1, drdt, dUdt_tmp, scratch, flag, vfrac);
            }

            if (m_advect_tracer) {
                compute_convective_fluxes_eb(lev, gbx, m_ntrac,
                                             fx, fy, fz, rhotrac, umac, vmac, wmac,
                                             get_tracer_bcrec().data(),
                                             get_tracer_bcrec_device_ptr(),
                                             flag, fcx, fcy, fcz, ccc);
                compute_convective_rate_eb(lev, gbx, m_ntrac, dUdt_tmp, fx, fy, fz,
                                           flag, vfrac, apx, apy, apz);
                redistribute_eb(lev, bx, m_ntrac, dtdt, dUdt_tmp, scratch, flag, vfrac);
            }
        }
        else
#endif
        {
            // velocity
            compute_convective_fluxes(lev, bx, AMREX_SPACEDIM, fx, fy, fz, vel,
                                      umac, vmac, wmac,
                                      get_velocity_bcrec().data(),
                                      get_velocity_bcrec_device_ptr());
            compute_convective_rate(lev, bx, AMREX_SPACEDIM, dvdt, fx, fy, fz);

            // density
            if (!m_constant_density) {
                compute_convective_fluxes(lev, bx, 1, fx, fy, fz, rho,
                                          umac, vmac, wmac,
                                          get_density_bcrec().data(),
                                          get_density_bcrec_device_ptr());
                compute_convective_rate(lev, bx, 1, drdt, fx, fy, fz);
            }

            // tracer
            if (m_advect_tracer) {
                compute_convective_fluxes(lev, bx, m_ntrac, fx, fy, fz, rhotrac,
                                          umac, vmac, wmac,
                                          get_tracer_bcrec().data(),
                                          get_tracer_bcrec_device_ptr());
                compute_convective_rate(lev, bx, m_ntrac, dtdt, fx, fy, fz);
            }
        }
    }
}
