#include <incflo.H>

using namespace amrex;

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
    int ngmac = nghost_mac();

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
            u_mac[lev].setBndry(0.0);
            v_mac[lev].setBndry(0.0);
            w_mac[lev].setBndry(0.0);
        }
        // Predict normal velocity to faces -- note that the {u_mac, v_mac, w_mac}
        //    returned from this call are on face CENTROIDS
        if (m_use_godunov) {
            predict_godunov(lev, time, u_mac[lev], v_mac[lev], w_mac[lev], *vel[lev], *density[lev],
                            *vel_forces[lev]);
        } else {
            predict_vels_on_faces(lev, time, u_mac[lev], v_mac[lev], w_mac[lev], *vel[lev]);
        }
    }

    apply_MAC_projection(u_mac, v_mac, w_mac, density, time);

    if (m_use_godunov) {
        VisMF::Write(u_mac[0], "u_mac");
        VisMF::Write(v_mac[0], "v_mac");
        VisMF::Write(w_mac[0], "w_mac");
        amrex::Abort("after predict mac and mac project");
    }

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
                                    (m_ntrac>0) ? conv_t[lev]->array(mfi) : Array4<Real>{},
                                    vel[lev]->const_array(mfi),
                                    density[lev]->array(mfi),
                                    (m_ntrac>0) ? tracer[lev]->array(mfi) : Array4<Real const>{},
                                    u_mac[lev].const_array(mfi),
                                    v_mac[lev].const_array(mfi),
                                    w_mac[lev].const_array(mfi));
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
                                 Array4<Real const> const& wmac)
{
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
        if (m_advect_tracer) {
            amrex::ParallelFor(bx, m_ntrac, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                dtdt(i,j,k,n) = 0.0;
            });
        }
        return;
    }

    bool regular = (flagfab.getType(amrex::grow(bx,1)) == FabType::regular);

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
    if (m_advect_tracer) nmaxcomp = std::max(nmaxcomp,m_ntrac);
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
    if (m_advect_tracer) {
        rhotracfab.resize(rhotrac_box, m_ntrac);
        eli_rt = rhotracfab.elixir();
        rhotrac = rhotracfab.array();
        amrex::ParallelFor(rhotrac_box, m_ntrac,
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
        if (!m_constant_density) {
            compute_convective_fluxes_eb(lev, gbx, 1, fx, fy, fz, rho, umac, vmac, wmac,
                                         get_density_bcrec().data(),
                                         get_density_bcrec_device_ptr(),
                                         flag, fcx, fcy, fcz, qface);
            compute_convective_rate_eb(lev, gbx, 1, dUdt_tmp, fx, fy, fz,
                                       flag, vfrac, apx, apy, apz);
            redistribute_eb(lev, bx, 1, drdt, dUdt_tmp, scratch, flag, vfrac);
        }

        if (m_advect_tracer) {
            compute_convective_fluxes_eb(lev, gbx, m_ntrac, fx, fy, fz, rhotrac, umac, vmac, wmac,
                                         get_tracer_bcrec().data(),
                                         get_tracer_bcrec_device_ptr(),
                                         flag, fcx, fcy, fcz, qface);
            compute_convective_rate_eb(lev, gbx, m_ntrac, dUdt_tmp, fx, fy, fz,
                                       flag, vfrac, apx, apy, apz);
            redistribute_eb(lev, bx, m_ntrac, dtdt, dUdt_tmp, scratch, flag, vfrac);
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
        if (!m_constant_density) {
            compute_convective_fluxes(lev, bx, 1, fx, fy, fz, rho, umac, vmac, wmac,
                                      get_density_bcrec().data(),
                                      get_density_bcrec_device_ptr());
            compute_convective_rate(lev, bx, 1, drdt, fx, fy, fz);
        }

        // tracer
        if (m_advect_tracer) {
            compute_convective_fluxes(lev, bx, m_ntrac, fx, fy, fz, rhotrac, umac, vmac, wmac,
                                      get_tracer_bcrec().data(),
                                      get_tracer_bcrec_device_ptr());
            compute_convective_rate(lev, bx, m_ntrac, dtdt, fx, fy, fz);
        }
    }
}
