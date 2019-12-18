#include <incflo_convection_K.H>

using namespace amrex;

void incflo_predict_vels_on_faces (Box const& ubx, Box const& vbx, Box const& wbx,
                                   Array4<Real> const& u, Array4<Real> const& v,
                                   Array4<Real> const& w, Array4<Real const> const& vcc)
{
    constexpr Real small = 1.e-10;

    amrex::ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real upls = vcc(i  ,j,k,0) - 0.5 * incflo_xslope(i  ,j,k,0,vcc);
        Real umns = vcc(i-1,j,k,0) + 0.5 * incflo_xslope(i-1,j,k,0,vcc);
        if (umns < 0.0 and upls > 0.0) {
            u(i,j,k) = 0.0;
        } else {
            Real avg = 0.5 * (upls + umns);
            if (std::abs(avg) < small) {
                u(i,j,k) = 0.0;
            } else if (avg > 0.0) {
                u(i,j,k) = umns;
            } else {
                u(i,j,k) = upls;
            }
        }
    });

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real vpls = vcc(i,j  ,k,1) - 0.5 * incflo_yslope(i,j  ,k,1,vcc);
        Real vmns = vcc(i,j-1,k,1) + 0.5 * incflo_yslope(i,j-1,k,1,vcc);
        if (vmns < 0.0 and vpls > 0.0) {
            v(i,j,k) = 0.0;
        } else {
            Real avg = 0.5 * (vpls + vmns);
            if (std::abs(avg) < small) {
                v(i,j,k) = 0.0;
            } else if (avg > 0.0) {
                v(i,j,k) = vmns;
            } else {
                v(i,j,k) = vpls;
            }
        }
    });

    amrex::ParallelFor(wbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real wpls = vcc(i,j,k  ,2) - 0.5 * incflo_zslope(i,j,k  ,2,vcc);
        Real wmns = vcc(i,j,k-1,2) + 0.5 * incflo_zslope(i,j,k-1,2,vcc);
        if (wmns < 0.0 and wpls > 0.0) {
            w(i,j,k) = 0.0;
        } else {
            Real avg = 0.5 * (wpls + wmns);
            if (std::abs(avg) < small) {
                w(i,j,k) = 0.0;
            } else if (avg > 0.0) {
                w(i,j,k) = wmns;
            } else {
                w(i,j,k) = wpls;
            }
        }
    });
}


#ifdef AMREX_USE_EB
void incflo_predict_vels_on_faces_eb (Box const& ccbx,
                                      Box const& ubx, Box const& vbx, Box const& wbx,
                                      Array4<Real> const& u, Array4<Real> const& v,
                                      Array4<Real> const& w, Array4<Real const> const& vcc,
                                      Array4<EBCellFlag const> const& flag,
                                      Array4<Real const> const& fcx,
                                      Array4<Real const> const& fcy,
                                      Array4<Real const> const& fcz)
{
    constexpr Real small = 1.e-10;

    FArrayBox tmp(amrex::grow(ccbx,1),2);
    Elixir eli = tmp.elixir();

    Array4<Real> const& upls = tmp.array(0);
    Array4<Real> const& umns = tmp.array(1);

    amrex::ParallelFor(amrex::grow(amrex::grow(ubx,1,1),2,1),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        upls(i,j,k) = vcc(i  ,j,k,0) - 0.5 * incflo_xslope(i  ,j,k,0,vcc);
        umns(i,j,k) = vcc(i-1,j,k,0) + 0.5 * incflo_xslope(i-1,j,k,0,vcc);
    });

    amrex::ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (flag(i,j,k).isConnected(-1,0,0)) {
            Real upls_on_centroid, umns_on_centroid;
            if (flag(i,j,k).isRegular()) {
                upls_on_centroid = upls(i,j,k);
                umns_on_centroid = umns(i,j,k);
            } else {
                int jj = j + static_cast<int>(std::copysign(1.0, fcx(i,j,k,0)));
                int kk = k + static_cast<int>(std::copysign(1.0, fcx(i,j,k,1)));

                Real fracy = std::abs(fcx(i,j,k,0));
                Real fracz = std::abs(fcx(i,j,k,1));

                upls_on_centroid = (1.0-fracy)*(1.0-fracz)*upls(i, j,k )+
                                        fracy *(1.0-fracz)*upls(i,jj,k )+
                                        fracz *(1.0-fracy)*upls(i, j,kk)+
                                        fracy *     fracz *upls(i,jj,kk);
                umns_on_centroid = (1.0-fracy)*(1.0-fracz)*umns(i, j,k )+
                                        fracy *(1.0-fracz)*umns(i,jj,k )+
                                        fracz *(1.0-fracy)*umns(i, j,kk)+
                                        fracy *     fracz *umns(i,jj,kk);
            }

            if (umns_on_centroid < 0.0 and upls_on_centroid > 0.0) {
                u(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (upls_on_centroid + umns_on_centroid);
                if (std::abs(avg) < small) {
                    u(i,j,k) = 0.0;
                } else if (avg > 0.0) {
                    u(i,j,k) = umns_on_centroid;
                } else {
                    u(i,j,k) = upls_on_centroid;
                }
            }
        } else {
            u(i,j,k) = 0.0;
        }
    });

    Array4<Real> const& vpls = upls;
    Array4<Real> const& vmns = umns;

    amrex::ParallelFor(amrex::grow(amrex::grow(vbx,0,1),2,1),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        vpls(i,j,k) = vcc(i,j  ,k,1) - 0.5 * incflo_yslope(i,j  ,k,1,vcc);
        vmns(i,j,k) = vcc(i,j-1,k,1) + 0.5 * incflo_yslope(i,j-1,k,1,vcc);
    });

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (flag(i,j,k).isConnected(0,-1,0)) {
            Real vpls_on_centroid, vmns_on_centroid;
            if (flag(i,j,k).isRegular()) {
                vpls_on_centroid = vpls(i,j,k);
                vmns_on_centroid = vmns(i,j,k);
            } else {
                int ii = i + static_cast<int>(std::copysign(1.0,fcy(i,j,k,0)));
                int kk = k + static_cast<int>(std::copysign(1.0,fcy(i,j,k,1)));

                Real fracx = std::abs(fcy(i,j,k,0));
                Real fracz = std::abs(fcy(i,j,k,1));

                vpls_on_centroid = (1.0-fracx)*(1.0-fracz)*vpls(i ,j,k )+
                                        fracx *(1.0-fracz)*vpls(ii,j,k )+
                                        fracz *(1.0-fracx)*vpls(i ,j,kk)+
                                        fracx *     fracz *vpls(ii,j,kk);
                vmns_on_centroid = (1.0-fracx)*(1.0-fracz)*vmns(i ,j,k )+
                                        fracx *(1.0-fracz)*vmns(ii,j,k )+
                                        fracz *(1.0-fracx)*vmns(i ,j,kk)+
                                        fracx *     fracz *vmns(ii,j,kk);
            }

            if (vmns_on_centroid < 0.0 and vpls_on_centroid > 0.0) {
                v(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (vpls_on_centroid + vmns_on_centroid);
                if (std::abs(avg) < small) {
                    v(i,j,k) = 0.0;
                } else if (avg > 0.0) {
                    v(i,j,k) = vmns_on_centroid;
                } else {
                    v(i,j,k) = vpls_on_centroid;
                }
            }
        } else {
            v(i,j,k) = 0.0;
        }
    });

    Array4<Real> const& wpls = upls;
    Array4<Real> const& wmns = umns;

    amrex::ParallelFor(amrex::grow(amrex::grow(vbx,0,1),1,1),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        wpls(i,j,k) = vcc(i,j,k  ,2) - 0.5 * incflo_zslope(i,j,k  ,2,vcc);
        wmns(i,j,k) = vcc(i,j,k-1,2) + 0.5 * incflo_zslope(i,j,k-1,2,vcc);
    });

    amrex::ParallelFor(wbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (flag(i,j,k).isConnected(0,0,-1)) {
            Real wpls_on_centroid, wmns_on_centroid;
            if (flag(i,j,k).isRegular()) {
                wpls_on_centroid = wpls(i,j,k);
                wmns_on_centroid = wmns(i,j,k);
            } else {
                int ii = i + static_cast<int>(std::copysign(1.0,fcz(i,j,k,0)));
                int jj = j + static_cast<int>(std::copysign(1.0,fcz(i,j,k,1)));

                Real fracx = std::abs(fcz(i,j,k,0));
                Real fracy = std::abs(fcz(i,j,k,1));

                wpls_on_centroid = (1.0-fracx)*(1.0-fracy)*wpls(i ,j ,k)+
                                        fracx *(1.0-fracy)*wpls(ii,j ,k)+
                                        fracy *(1.0-fracx)*wpls(i ,jj,k)+
                                        fracx *     fracy *wpls(ii,jj,k);
                wmns_on_centroid = (1.0-fracx)*(1.0-fracy)*wmns(i ,j ,k)+
                                        fracx *(1.0-fracy)*wmns(ii,j ,k)+
                                        fracy *(1.0-fracx)*wmns(i ,jj,k)+
                                        fracx *     fracy *wmns(ii,jj,k);
            }

            if (wmns_on_centroid < 0.0 and wpls_on_centroid > 0.0) {
                w(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (wpls_on_centroid + wmns_on_centroid);
                if (std::abs(avg) < small) {
                    w(i,j,k) = 0.0;
                } else if (avg > 0.0) {
                    w(i,j,k) = wmns_on_centroid;
                } else {
                    w(i,j,k) = wpls_on_centroid;
                }
            }
        } else {
            w(i,j,k) = 0.0;
        }
    });
}
#endif
