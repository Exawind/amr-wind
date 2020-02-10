#include <incflo.H>
#include <incflo_MAC_bcs.H>

using namespace amrex;

void incflo::predict_vels_on_faces (int lev, Real time, MultiFab& u_mac, MultiFab& v_mac,
                                    MultiFab& w_mac, MultiFab const& vel)
{
#ifdef AMREX_USE_EB
    auto const& fact = this->EBFactory(lev);
    auto const& flags = fact.getMultiEBCellFlagFab();
    auto const& fcent = fact.getFaceCent();
    auto const& ccent = fact.getCentroid();
#endif

    Box const& domain = Geom(lev).Domain();
    Vector<BCRec> const& h_bcrec = get_velocity_bcrec();
    BCRec const* d_bcrec = get_velocity_bcrec_device_ptr();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(vel, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& ubx = mfi.nodaltilebox(0);
            Box const& vbx = mfi.nodaltilebox(1);
            Box const& wbx = mfi.nodaltilebox(2);
            Array4<Real> const& u = u_mac.array(mfi);
            Array4<Real> const& v = v_mac.array(mfi);
            Array4<Real> const& w = w_mac.array(mfi);
            Array4<Real const> const& vcc = vel.const_array(mfi);
#ifdef AMREX_USE_EB
            Box const& bx = mfi.tilebox();
            EBCellFlagFab const& flagfab = flags[mfi];
            Array4<EBCellFlag const> const& flagarr = flagfab.const_array();
            auto const typ = flagfab.getType(amrex::grow(bx,1));
            if (typ == FabType::covered)
            {
                amrex::ParallelFor(ubx, vbx, wbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { u(i,j,k) = 0.0; },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { v(i,j,k) = 0.0; },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { w(i,j,k) = 0.0; });
            }
            else if (typ == FabType::singlevalued)
            {
                Array4<Real const> const& fcx = fcent[0]->const_array(mfi);
                Array4<Real const> const& fcy = fcent[1]->const_array(mfi);
                Array4<Real const> const& fcz = fcent[2]->const_array(mfi);
                Array4<Real const> const& ccc = ccent.const_array(mfi);
                predict_vels_on_faces_eb(lev,bx,ubx,vbx,wbx,u,v,w,vcc,flagarr,fcx,fcy,fcz,ccc);
            }
            else
#endif
            {
                predict_vels_on_faces(lev,ubx,vbx,wbx,u,v,w,vcc);
            }

            incflo_set_mac_bcs(domain,ubx,vbx,wbx,u,v,w,vcc,h_bcrec,d_bcrec);
        }
    }
}
