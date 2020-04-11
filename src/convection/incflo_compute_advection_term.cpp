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
    int ngmac = nghost_mac();

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Predict normal velocity to faces -- note that the {u_mac, v_mac, w_mac}
        //    returned from this call are on face CENTROIDS
        if (m_use_godunov) {
            predict_godunov(lev, time, *u_mac[lev], *v_mac[lev], *w_mac[lev], *vel[lev], *vel_forces[lev]);
        } else {
            predict_vels_on_faces(lev, *u_mac[lev], *v_mac[lev], *w_mac[lev], *vel[lev]);
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
incflo::compute_convective_term (Box const& bx, int lev,
                                 MFIter const&,
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

    Box rhotrac_box = amrex::grow(bx,2);
    if (m_use_godunov) rhotrac_box.grow(1);

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
                                  velocity().bcrec_device().data(),
                                  get_velocity_iconserv_device_ptr(),
                                  tmpfab.dataPtr());
        if (!m_constant_density) {
            compute_godunov_advection(lev, bx, 1,
                                      drdt, rho,
                                      umac, vmac, wmac, {},
                                      density().bcrec_device().data(),
                                      get_density_iconserv_device_ptr(),
                                      tmpfab.dataPtr());
        }
        if (m_advect_tracer) {
            compute_godunov_advection(lev, bx, m_ntrac,
                                      dtdt, rhotrac,
                                      umac, vmac, wmac, ftra,
                                      tracer().bcrec_device().data(),
                                      get_tracer_iconserv_device_ptr(),
                                      tmpfab.dataPtr());
        }
        Gpu::streamSynchronize();
    }
    else
    {
        Box tmpbox = amrex::surroundingNodes(bx);
        int tmpcomp = nmaxcomp*AMREX_SPACEDIM;

        FArrayBox tmpfab(tmpbox, tmpcomp);
        Elixir eli = tmpfab.elixir();

        Array4<Real> fx = tmpfab.array(0);
        Array4<Real> fy = tmpfab.array(nmaxcomp);
        Array4<Real> fz = tmpfab.array(nmaxcomp*2);

        // velocity
        compute_convective_fluxes(lev, bx, AMREX_SPACEDIM, fx, fy, fz, vel,
                                  umac, vmac, wmac,
                                  velocity().bcrec().data(),
                                  velocity().bcrec_device().data());
        compute_convective_rate(lev, bx, AMREX_SPACEDIM, dvdt, fx, fy, fz);

        // density
        if (!m_constant_density) {
            compute_convective_fluxes(lev, bx, 1, fx, fy, fz, rho,
                                      umac, vmac, wmac,
                                      density().bcrec().data(),
                                      density().bcrec_device().data());
            compute_convective_rate(lev, bx, 1, drdt, fx, fy, fz);
        }

        // tracer
        if (m_advect_tracer) {
            compute_convective_fluxes(lev, bx, m_ntrac, fx, fy, fz, rhotrac,
                                      umac, vmac, wmac,
                                      tracer().bcrec().data(),
                                      tracer().bcrec_device().data());
            compute_convective_rate(lev, bx, m_ntrac, dtdt, fx, fy, fz);
        }
        
    }
}

void incflo::compute_convective_rate (int lev, Box const& bx, int ncomp,
                                      Array4<Real> const& dUdt,
                                      Array4<Real const> const& fx,
                                      Array4<Real const> const& fy,
                                      Array4<Real const> const& fz)
{
    const auto dxinv = Geom(lev).InvCellSizeArray();
    amrex::ParallelFor(bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        dUdt(i,j,k,n) = dxinv[0] * (fx(i,j,k,n) - fx(i+1,j,k,n))
            +           dxinv[1] * (fy(i,j,k,n) - fy(i,j+1,k,n))
            +           dxinv[2] * (fz(i,j,k,n) - fz(i,j,k+1,n));
    });
}
