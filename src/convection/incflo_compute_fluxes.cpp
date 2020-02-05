#include <incflo_convection_K.H>
#include <incflo.H>
#include <utility>

using namespace amrex;

namespace {
    std::pair<bool,bool> has_extdir (BCRec const* bcrec, int ncomp, int dir)
    {
        std::pair<bool,bool> r{false,false};
        for (int n = 0; n < ncomp; ++n) {
            r.first = r.first or bcrec[n].lo(dir) == BCType::ext_dir;
            r.second = r.second or bcrec[n].hi(dir) == BCType::ext_dir;
        }
        return r;
    }
}

void
incflo::compute_convective_fluxes (int lev, Box const& bx, int ncomp,
                                   Array4<Real> const& fx,
                                   Array4<Real> const& fy,
                                   Array4<Real> const& fz,
                                   Array4<Real const> const& q,
                                   Array4<Real const> const& umac,
                                   Array4<Real const> const& vmac,
                                   Array4<Real const> const& wmac,
                                   BCRec const* h_bcrec, BCRec const* d_bcrec)
{
    constexpr Real small = 1.e-10;

    const Box& domain_box = geom[lev].Domain();
    const int domain_ilo = domain_box.smallEnd(0);
    const int domain_ihi = domain_box.bigEnd(0);
    const int domain_jlo = domain_box.smallEnd(1);
    const int domain_jhi = domain_box.bigEnd(1);
    const int domain_klo = domain_box.smallEnd(2);
    const int domain_khi = domain_box.bigEnd(2);

    Box const& xbx = amrex::surroundingNodes(bx,0);
    Box const& ybx = amrex::surroundingNodes(bx,1);
    Box const& zbx = amrex::surroundingNodes(bx,2);

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    auto extdir_lohi = has_extdir(h_bcrec, ncomp, static_cast<int>(Direction::x));
    bool has_extdir_lo = extdir_lohi.first;
    bool has_extdir_hi = extdir_lohi.second;
    if ((has_extdir_lo and domain_ilo >= xbx.smallEnd(0)-1) or
        (has_extdir_hi and domain_ihi <= xbx.bigEnd(0)))
    {
        amrex::ParallelFor(xbx, ncomp, [d_bcrec,q,domain_ilo,domain_ihi,umac,small,fx]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            bool extdir_ilo = d_bcrec[n].lo(0) == BCType::ext_dir;
            bool extdir_ihi = d_bcrec[n].hi(0) == BCType::ext_dir;
            Real qs;
            if (extdir_ilo and i <= domain_ilo) {
                qs = q(domain_ilo-1,j,k,n);
            } else if (extdir_ihi and i >= domain_ihi+1) {
                qs = q(domain_ihi+1,j,k,n);
            } else {
                Real qpls = q(i,j,k,n) - 0.5 * incflo_xslope_extdir
                    (i,j,k,n,q, extdir_ilo, extdir_ihi, domain_ilo, domain_ihi);
                Real qmns = q(i-1,j,k,n) + 0.5 * incflo_xslope_extdir
                    (i-1,j,k,n,q, extdir_ilo, extdir_ihi, domain_ilo, domain_ihi);
                if (umac(i,j,k) > small) {
                    qs = qmns;
                } else if (umac(i,j,k) < -small) {
                    qs = qpls;
                } else {
                    qs = 0.5*(qmns+qpls);
                }
            }
            fx(i,j,k,n) = qs * umac(i,j,k);
        });
    }
    else
    {
        amrex::ParallelFor(xbx, ncomp, [q,umac,small,fx]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qpls = q(i  ,j,k,n) - 0.5 * incflo_xslope(i  ,j,k,n,q);
            Real qmns = q(i-1,j,k,n) + 0.5 * incflo_xslope(i-1,j,k,n,q);
            Real qs;
            if (umac(i,j,k) > small) {
                qs = qmns;
            } else if (umac(i,j,k) < -small) {
                qs = qpls;
            } else {
                qs = 0.5*(qmns+qpls);
            }
            fx(i,j,k,n) = qs * umac(i,j,k);
        });
    }

    extdir_lohi = has_extdir(h_bcrec, ncomp,  static_cast<int>(Direction::y));
    has_extdir_lo = extdir_lohi.first;
    has_extdir_hi = extdir_lohi.second;
    if ((has_extdir_lo and domain_jlo >= ybx.smallEnd(1)-1) or
        (has_extdir_hi and domain_jhi <= ybx.bigEnd(1)))
    {
        amrex::ParallelFor(ybx, ncomp, [d_bcrec,q,domain_jlo,domain_jhi,vmac,small,fy]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            bool extdir_jlo = d_bcrec[n].lo(1) == BCType::ext_dir;
            bool extdir_jhi = d_bcrec[n].hi(1) == BCType::ext_dir;
            Real qs;
            if (extdir_jlo and j <= domain_jlo) {
                qs = q(i,domain_jlo-1,k,n);
            } else if (extdir_jhi and j >= domain_jhi+1) {
                qs = q(i,domain_jhi+1,k,n);
            } else {
                Real qpls = q(i,j,k,n) - 0.5 * incflo_yslope_extdir
                    (i,j,k,n,q, extdir_jlo, extdir_jhi, domain_jlo, domain_jhi);
                Real qmns = q(i,j-1,k,n) + 0.5 * incflo_yslope_extdir
                    (i,j-1,k,n,q, extdir_jlo, extdir_jhi, domain_jlo, domain_jhi);
                if (vmac(i,j,k) > small) {
                    qs = qmns;
                } else if (vmac(i,j,k) < -small) {
                    qs = qpls;
                } else {
                    qs = 0.5*(qmns+qpls);
                }
            }
            fy(i,j,k,n) = qs * vmac(i,j,k);
        });
    }
    else
    {
        amrex::ParallelFor(ybx, ncomp, [q,vmac,small,fy]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qpls = q(i,j  ,k,n) - 0.5 * incflo_yslope(i,j  ,k,n,q);
            Real qmns = q(i,j-1,k,n) + 0.5 * incflo_yslope(i,j-1,k,n,q);
            Real qs;
            if (vmac(i,j,k) > small) {
                qs = qmns;
            } else if (vmac(i,j,k) < -small) {
                qs = qpls;
            } else {
                qs = 0.5*(qmns+qpls);
            }
            fy(i,j,k,n) = qs * vmac(i,j,k);
        });
    }

    extdir_lohi = has_extdir(h_bcrec, ncomp, static_cast<int>(Direction::z));
    has_extdir_lo = extdir_lohi.first;
    has_extdir_hi = extdir_lohi.second;
    if ((has_extdir_lo and domain_klo >= zbx.smallEnd(2)-1) or
        (has_extdir_hi and domain_khi <= zbx.bigEnd(2)))
    {
        amrex::ParallelFor(zbx, ncomp, [d_bcrec,q,domain_klo,domain_khi,wmac,small,fz]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            bool extdir_klo = d_bcrec[n].lo(2) == BCType::ext_dir;
            bool extdir_khi = d_bcrec[n].hi(2) == BCType::ext_dir;
            Real qs;
            if (extdir_klo and k <= domain_klo) {
                qs = q(i,j,domain_klo-1,n);
            } else if (extdir_khi and k >= domain_khi+1) {
                qs = q(i,j,domain_khi+1,n);
            } else {
                Real qpls = q(i,j,k,n) - 0.5 * incflo_zslope_extdir
                    (i,j,k,n,q, extdir_klo, extdir_khi, domain_klo, domain_khi);
                Real qmns = q(i,j,k-1,n) + 0.5 * incflo_zslope_extdir(
                    i,j,k-1,n,q, extdir_klo, extdir_khi, domain_klo, domain_khi);
                if (wmac(i,j,k) > small) {
                    qs = qmns;
                } else if (wmac(i,j,k) < -small) {
                    qs = qpls;
                } else {
                    qs = 0.5*(qmns+qpls);
                }
            }
            fz(i,j,k,n) = qs * wmac(i,j,k);
        });
    }
    else
    {
        amrex::ParallelFor(zbx, ncomp, [q,wmac,small,fz]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qpls = q(i,j,k  ,n) - 0.5 * incflo_zslope(i,j,k  ,n,q);
            Real qmns = q(i,j,k-1,n) + 0.5 * incflo_zslope(i,j,k-1,n,q);
            Real qs;
            if (wmac(i,j,k) > small) {
                qs = qmns;
            } else if (wmac(i,j,k) < -small) {
                qs = qpls;
            } else {
                qs = 0.5*(qmns+qpls);
            }
            fz(i,j,k,n) = qs * wmac(i,j,k);
        });
    }
}

#ifdef AMREX_USE_EB
void incflo::compute_convective_fluxes_eb (int lev, Box const& bx, int ncomp,
                                           Array4<Real> const& fx,
                                           Array4<Real> const& fy,
                                           Array4<Real> const& fz,
                                           Array4<Real const> const& q,
                                           Array4<Real const> const& umac,
                                           Array4<Real const> const& vmac,
                                           Array4<Real const> const& wmac,
                                           BCRec const* h_bcrec,
                                           BCRec const* d_bcrec,
                                           Array4<EBCellFlag const> const& flag,
                                           Array4<Real const> const& fcx,
                                           Array4<Real const> const& fcy,
                                           Array4<Real const> const& fcz,
                                           Array4<Real> const& qface)
{
    constexpr Real small = 1.e-10;

    const Box& domain_box = geom[lev].Domain();
    const int domain_ilo = domain_box.smallEnd(0);
    const int domain_ihi = domain_box.bigEnd(0);
    const int domain_jlo = domain_box.smallEnd(1);
    const int domain_jhi = domain_box.bigEnd(1);
    const int domain_klo = domain_box.smallEnd(2);
    const int domain_khi = domain_box.bigEnd(2);

    Box const& xbx = amrex::surroundingNodes(bx,0);
    Box const& ybx = amrex::surroundingNodes(bx,1);
    Box const& zbx = amrex::surroundingNodes(bx,2);

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    auto extdir_lohi = has_extdir(h_bcrec, ncomp, static_cast<int>(Direction::x));
    bool has_extdir_lo = extdir_lohi.first;
    bool has_extdir_hi = extdir_lohi.second;
    if ((has_extdir_lo and domain_ilo >= xbx.smallEnd(0)-1) or
        (has_extdir_hi and domain_ihi <= xbx.bigEnd(0)))
    {
        amrex::ParallelFor(Box(xbx).grow(1,1).grow(2,1), ncomp,
        [d_bcrec,domain_ilo,domain_ihi,q,flag,umac,small,qface] 
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            bool extdir_ilo = d_bcrec[n].lo(0) == BCType::ext_dir;
            bool extdir_ihi = d_bcrec[n].hi(0) == BCType::ext_dir;
            Real qs;
            if (extdir_ilo and i <= domain_ilo) {
                qs = q(domain_ilo-1,j,k,n);
            } else if (extdir_ihi and i >= domain_ihi+1) {
                qs = q(domain_ihi+1,j,k,n);
            } else {
                Real qpls = q(i,j,k,n) - 0.5 * incflo_xslope_extdir_eb
                    (i,j,k,n,q,flag, extdir_ilo, extdir_ihi, domain_ilo, domain_ihi);
                Real qmns = q(i-1,j,k,n) + 0.5 * incflo_xslope_extdir_eb
                    (i-1,j,k,n,q,flag, extdir_ilo, extdir_ihi, domain_ilo, domain_ihi);
                if (umac(i,j,k) > small) {
                    qs = qmns;
                } else if (umac(i,j,k) < -small) {
                    qs = qpls;
                } else {
                    qs = 0.5*(qmns+qpls);
                }
            }
            qface(i,j,k,n) = qs;
        });
    }
    else
    {
        amrex::ParallelFor(Box(xbx).grow(1,1).grow(2,1), ncomp,
        [q,flag,umac,small,qface]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qpls = q(i  ,j,k,n) - 0.5 * incflo_xslope_eb(i  ,j,k,n,q,flag);
            Real qmns = q(i-1,j,k,n) + 0.5 * incflo_xslope_eb(i-1,j,k,n,q,flag);
            Real qs;
            if (umac(i,j,k) > small) {
                qs = qmns;
            } else if (umac(i,j,k) < -small) {
                qs = qpls;
            } else {
                qs = 0.5*(qmns+qpls);
            }
            qface(i,j,k,n) = qs;
        });
    }

    amrex::ParallelFor(xbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isConnected(-1,0,0)) {
            Real qcent;
            if (flag(i,j,k).isRegular() or flag(i-1,j,k).isRegular()) {
                qcent = qface(i,j,k,n);
            } else {
                int jj = j + static_cast<int>(std::copysign(1.0, fcx(i,j,k,0)));
                int kk = k + static_cast<int>(std::copysign(1.0, fcx(i,j,k,1)));

                Real fracy = std::abs(fcx(i,j,k,0));
                Real fracz = std::abs(fcx(i,j,k,1));

                qcent = (1.0-fracy)*(1.0-fracz)*qface(i, j,k ,n)+
                             fracy *(1.0-fracz)*qface(i,jj,k ,n)+
                             fracz *(1.0-fracy)*qface(i, j,kk,n)+
                             fracy *     fracz *qface(i,jj,kk,n);
            }
            fx(i,j,k,n) = qcent * umac(i,j,k);
        } else {
            fx(i,j,k,n) = 0.0;
        }
    });

    extdir_lohi = has_extdir(h_bcrec, ncomp,  static_cast<int>(Direction::y));
    has_extdir_lo = extdir_lohi.first;
    has_extdir_hi = extdir_lohi.second;
    if ((has_extdir_lo and domain_jlo >= ybx.smallEnd(1)-1) or
        (has_extdir_hi and domain_jhi <= ybx.bigEnd(1)))
    {
        amrex::ParallelFor(Box(ybx).grow(0,1).grow(2,1), ncomp,
        [d_bcrec,domain_jlo,domain_jhi,q,flag,vmac,small,qface]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            bool extdir_jlo = d_bcrec[n].lo(1) == BCType::ext_dir;
            bool extdir_jhi = d_bcrec[n].hi(1) == BCType::ext_dir;
            Real qs;
            if (extdir_jlo and j <= domain_jlo) {
                qs = q(i,domain_jlo-1,k,n);
            } else if (extdir_jhi and j >= domain_jhi+1) {
                qs = q(i,domain_jhi+1,k,n);
            } else {
                Real qpls = q(i,j,k,n) - 0.5 * incflo_yslope_extdir_eb
                    (i,j,k,n,q,flag, extdir_jlo, extdir_jhi, domain_jlo, domain_jhi);
                Real qmns = q(i,j-1,k,n) + 0.5 * incflo_yslope_extdir_eb
                    (i,j-1,k,n,q,flag, extdir_jlo, extdir_jhi, domain_jlo, domain_jhi);
                if (vmac(i,j,k) > small) {
                    qs = qmns;
                } else if (vmac(i,j,k) < -small) {
                    qs = qpls;
                } else {
                    qs = 0.5*(qmns+qpls);
                }
            }
            qface(i,j,k,n) = qs;
        });
    }
    else
    {
        amrex::ParallelFor(Box(ybx).grow(0,1).grow(2,1), ncomp,
        [q,flag,vmac,small,qface]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qpls = q(i,j  ,k,n) - 0.5 * incflo_yslope_eb(i,j  ,k,n,q,flag);
            Real qmns = q(i,j-1,k,n) + 0.5 * incflo_yslope_eb(i,j-1,k,n,q,flag);
            Real qs;
            if (vmac(i,j,k) > small) {
                qs = qmns;
            } else if (vmac(i,j,k) < -small) {
                qs = qpls;
            } else {
                qs = 0.5*(qmns+qpls);
            }
            qface(i,j,k,n) = qs;
        });
    }

    amrex::ParallelFor(ybx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isConnected(0,-1,0)) {
            Real qcent;
            if (flag(i,j,k).isRegular() or flag(i,j-1,k).isRegular()) {
                qcent = qface(i,j,k,n);
            } else {
                int ii = i + static_cast<int>(std::copysign(1.0,fcy(i,j,k,0)));
                int kk = k + static_cast<int>(std::copysign(1.0,fcy(i,j,k,1)));

                Real fracx = std::abs(fcy(i,j,k,0));
                Real fracz = std::abs(fcy(i,j,k,1));

                qcent = (1.0-fracx)*(1.0-fracz)*qface(i ,j,k ,n)+
                             fracx *(1.0-fracz)*qface(ii,j,k ,n)+
                             fracz *(1.0-fracx)*qface(i ,j,kk,n)+
                             fracx *     fracz *qface(ii,j,kk,n);

            }
            fy(i,j,k,n) = qcent * vmac(i,j,k);
        } else {
            fy(i,j,k,n) = 0.0;
        }
    });

    extdir_lohi = has_extdir(h_bcrec, ncomp, static_cast<int>(Direction::z));
    has_extdir_lo = extdir_lohi.first;
    has_extdir_hi = extdir_lohi.second;
    if ((has_extdir_lo and domain_klo >= zbx.smallEnd(2)-1) or
        (has_extdir_hi and domain_khi <= zbx.bigEnd(2)))
    {
        amrex::ParallelFor(Box(zbx).grow(0,1).grow(1,1), ncomp,
        [d_bcrec,domain_klo,domain_khi,q,flag,wmac,small,qface]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            bool extdir_klo = d_bcrec[n].lo(2) == BCType::ext_dir;
            bool extdir_khi = d_bcrec[n].hi(2) == BCType::ext_dir;
            Real qs;
            if (extdir_klo and k <= domain_klo) {
                qs = q(i,j,domain_klo-1,n);
            } else if (extdir_khi and k >= domain_khi+1) {
                qs = q(i,j,domain_khi+1,n);
            } else {
                Real qpls = q(i,j,k,n) - 0.5 * incflo_zslope_extdir_eb
                    (i,j,k,n,q,flag, extdir_klo, extdir_khi, domain_klo, domain_khi);
                Real qmns = q(i,j,k-1,n) + 0.5 * incflo_zslope_extdir_eb(
                    i,j,k-1,n,q,flag, extdir_klo, extdir_khi, domain_klo, domain_khi);
                if (wmac(i,j,k) > small) {
                    qs = qmns;
                } else if (wmac(i,j,k) < -small) {
                    qs = qpls;
                } else {
                    qs = 0.5*(qmns+qpls);
                }
            }
            qface(i,j,k,n) = qs;
        });
    }
    else
    {
        amrex::ParallelFor(Box(zbx).grow(0,1).grow(1,1), ncomp,
        [q,flag,wmac,small,qface]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qpls = q(i,j,k  ,n) - 0.5 * incflo_zslope_eb(i,j,k  ,n,q,flag);
            Real qmns = q(i,j,k-1,n) + 0.5 * incflo_zslope_eb(i,j,k-1,n,q,flag);
            Real qs;
            if (wmac(i,j,k) > small) {
                qs = qmns;
            } else if (wmac(i,j,k) < -small) {
                qs = qpls;
            } else {
                qs = 0.5*(qmns+qpls);
            }
            qface(i,j,k,n) = qs;
        });
    }

    amrex::ParallelFor(zbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isConnected(0,0,-1)) {
            Real qcent;
            if (flag(i,j,k).isRegular() or flag(i,j,k-1).isRegular()) {
                qcent = qface(i,j,k,n);
            } else {
                int ii = i + static_cast<int>(std::copysign(1.0,fcz(i,j,k,0)));
                int jj = j + static_cast<int>(std::copysign(1.0,fcz(i,j,k,1)));

                Real fracx = std::abs(fcz(i,j,k,0));
                Real fracy = std::abs(fcz(i,j,k,1));

                qcent = (1.0-fracx)*(1.0-fracy)*qface(i ,j ,k,n)+
                             fracx *(1.0-fracy)*qface(ii,j ,k,n)+
                             fracy *(1.0-fracx)*qface(i ,jj,k,n)+
                             fracx *     fracy *qface(ii,jj,k,n);
            }
            fz(i,j,k,n) = qcent * wmac(i,j,k);
        } else {
            fz(i,j,k,n) = 0.0;
        }
    });
}
#endif
