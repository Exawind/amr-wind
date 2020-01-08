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
        amrex::ParallelFor(xbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
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
        amrex::ParallelFor(xbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
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
        amrex::ParallelFor(ybx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
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
        amrex::ParallelFor(ybx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
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
        amrex::ParallelFor(zbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
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
        amrex::ParallelFor(zbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
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
    bool has_extdir_lo = (h_bcrec[0].lo(0) == BCType::ext_dir) or
                         (h_bcrec[1].lo(0) == BCType::ext_dir) or
                         (h_bcrec[2].lo(0) == BCType::ext_dir);
    bool has_extdir_hi = (h_bcrec[0].hi(0) == BCType::ext_dir) or
                         (h_bcrec[1].hi(0) == BCType::ext_dir) or
                         (h_bcrec[2].hi(0) == BCType::ext_dir);
    if ((has_extdir_lo and domain_ilo >= xbx.smallEnd(0)-1) or
        (has_extdir_hi and domain_ihi <= xbx.bigEnd(0)))
    {
        amrex::ParallelFor(amrex::grow(amrex::grow(xbx,1,1),2,1), ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
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
        amrex::ParallelFor(amrex::grow(amrex::grow(xbx,1,1),2,1), ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
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

    has_extdir_lo = (h_bcrec[0].lo(1) == BCType::ext_dir) or
                    (h_bcrec[1].lo(1) == BCType::ext_dir) or
                    (h_bcrec[2].lo(1) == BCType::ext_dir);
    has_extdir_hi = (h_bcrec[0].hi(1) == BCType::ext_dir) or
                    (h_bcrec[1].hi(1) == BCType::ext_dir) or
                    (h_bcrec[2].hi(1) == BCType::ext_dir);
    if ((has_extdir_lo and domain_jlo >= ybx.smallEnd(1)-1) or
        (has_extdir_hi and domain_jhi <= ybx.bigEnd(1)))
    {
        amrex::ParallelFor(amrex::grow(amrex::grow(ybx,0,1),2,1), ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
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
        amrex::ParallelFor(amrex::grow(amrex::grow(ybx,0,1),2,1), ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
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

    has_extdir_lo = (h_bcrec[0].lo(2) == BCType::ext_dir) or
                    (h_bcrec[1].lo(2) == BCType::ext_dir) or
                    (h_bcrec[2].lo(2) == BCType::ext_dir);
    has_extdir_hi = (h_bcrec[0].hi(2) == BCType::ext_dir) or
                    (h_bcrec[1].hi(2) == BCType::ext_dir) or
                    (h_bcrec[2].hi(2) == BCType::ext_dir);
    if ((has_extdir_lo and domain_klo >= zbx.smallEnd(2)-1) or
        (has_extdir_hi and domain_khi <= zbx.bigEnd(2)))
    {
        amrex::ParallelFor(amrex::grow(amrex::grow(zbx,0,1),1,1), ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
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
        amrex::ParallelFor(amrex::grow(amrex::grow(zbx,0,1),1,1), ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
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


namespace ugradu_aux {

//
// Compute upwind non-normal velocity
//
AMREX_GPU_HOST_DEVICE
Real
upwind(const Real velocity_minus,
       const Real velocity_plus,
       const Real u_edge)
{
  // Small value to protect against tiny velocities used in upwinding
  const Real small_velocity(1.e-10);

  if(std::abs(u_edge) < small_velocity)
    return .5*(velocity_minus+velocity_plus);

  return u_edge > 0 ? velocity_minus : velocity_plus;
}

AMREX_GPU_HOST_DEVICE
bool
is_equal_to_any(const int bc,
                const int* bc_types,
                const int size)
{
  for(int i(0); i < size; ++i)
  {
    if(bc == bc_types[i])
      return true;
  }
  return false;
}

} // end namespace ugradu_aux

using namespace ugradu_aux;

//
// Compute the three components of the convection term
//
void
incflo::incflo_compute_fluxes(int lev,
                          Vector< std::unique_ptr<MultiFab> >& a_fx,
                          Vector< std::unique_ptr<MultiFab> >& a_fy,
                          Vector< std::unique_ptr<MultiFab> >& a_fz,
                          Vector< std::unique_ptr<MultiFab> >& state_in,
                          const int state_comp, const int ncomp,
                          Vector< std::unique_ptr<MultiFab> >& xslopes_in,
                          Vector< std::unique_ptr<MultiFab> >& yslopes_in,
                          Vector< std::unique_ptr<MultiFab> >& zslopes_in,
                          const int slopes_comp,
                          Vector< std::unique_ptr<MultiFab> >& u_mac,
                          Vector< std::unique_ptr<MultiFab> >& v_mac,
                          Vector< std::unique_ptr<MultiFab> >& w_mac)
{
        Box domain(geom[lev].Domain());

#ifdef AMREX_USE_EB
        // Get EB geometric info
        Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
        Array< const MultiCutFab*,AMREX_SPACEDIM> facecent;
        const amrex::MultiFab*                    volfrac;
        const amrex::MultiCutFab*                 bndrycent;

        areafrac  =   ebfactory[lev] -> getAreaFrac();
        facecent  =   ebfactory[lev] -> getFaceCent();
        volfrac   = &(ebfactory[lev] -> getVolFrac());
        bndrycent = &(ebfactory[lev] -> getBndryCent());
#endif

        // Create cc_mask
        iMultiFab cc_mask(grids[lev], dmap[lev], 1, 1);

        const int covered_value = 1;
        const int notcovered_value = 0;
        const int physical_boundaries_value = 0;
        const int interior_value = 1;

        cc_mask.BuildMask(geom[lev].Domain(), geom[lev].periodicity(),
            covered_value, notcovered_value,
            physical_boundaries_value, interior_value);

#ifdef AMREX_USE_EB
        // We do this here to avoid any confusion about the FAB setVal.
        a_fx[lev]->setVal(covered_val, 0, ncomp);
        a_fy[lev]->setVal(covered_val, 0, ncomp);
        a_fz[lev]->setVal(covered_val, 0, ncomp);
#endif

        for (MFIter mfi(*state_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox ();

#ifdef AMREX_USE_EB
            // this is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox&  state_fab = static_cast<EBFArrayBox const&>((*state_in[lev])[mfi]);
            const EBCellFlagFab&  flags = state_fab.getEBCellFlagFab();

            if (flags.getType(amrex::grow(bx,0)) != FabType::covered )
            {
                // No cut cells in tile + 1 halo -> use non-eb routine
                if (flags.getType(amrex::grow(bx,1)) == FabType::regular )
                {
                    incflo_compute_fluxes_on_box(lev, bx, (*a_fx[lev])[mfi], (*a_fy[lev])[mfi], (*a_fz[lev])[mfi],
                                                 (*state_in[lev])[mfi], state_comp, ncomp,
                                                 (*xslopes_in[lev])[mfi], (*yslopes_in[lev])[mfi], (*zslopes_in[lev])[mfi], slopes_comp,
                                                 (*u_mac[lev])[mfi], (*v_mac[lev])[mfi], (*w_mac[lev])[mfi]);

                }
                else
                {
                    incflo_compute_eb_fluxes_on_box(lev, bx, (*a_fx[lev])[mfi], (*a_fy[lev])[mfi], (*a_fz[lev])[mfi],
                                                    (*state_in[lev])[mfi], state_comp, ncomp,
                                                    (*xslopes_in[lev])[mfi], (*yslopes_in[lev])[mfi], (*zslopes_in[lev])[mfi], slopes_comp,
                                                    ( *u_mac[lev])[mfi], ( *v_mac[lev])[mfi], ( *w_mac[lev])[mfi],
                                                    (*areafrac[0])[mfi], (*areafrac[1])[mfi], (*areafrac[2])[mfi],
                                                    (*facecent[0])[mfi], (*facecent[1])[mfi], (*facecent[2])[mfi],
                                                    (*volfrac)[mfi], (*bndrycent)[mfi], cc_mask[mfi], flags);
                      }
            }
#else
                    // HACK HACK HACK -- NOTE THIS IS NOT READY FOR USE
                    if (use_godunov)
                    {
                       // These are place-holders for now
                       MultiFab tforces(grids[lev], dmap[lev], ncomp, 2);;
                       MultiFab divu   (grids[lev], dmap[lev],     1, 2);;
                       tforces.setVal(0.0);
                       divu.setVal(0.0);

                       BCRec dom_bc;
                       {
                         // const int* lo_bc = phys_bc.lo();
                         // const int* hi_bc = phys_bc.hi();
                         // HACK -- just set to all int_dir as stand-in for periodic 
                         dom_bc.setLo(0,BCType::int_dir);
                         dom_bc.setHi(0,BCType::int_dir);
                         dom_bc.setLo(1,BCType::int_dir);
                         dom_bc.setHi(1,BCType::int_dir);
                         dom_bc.setLo(2,BCType::int_dir);
                         dom_bc.setHi(2,BCType::int_dir);
                       }

                       Gpu::ManagedVector<BCRec> bc(ncomp);
                       for (int n = 0; n < ncomp; ++n) 
                            setBC(bx, geom[lev].Domain(), dom_bc, bc[n]);

                       incflo_godunov_fluxes_on_box(lev, bx, (*a_fx[lev]).array(mfi), (*a_fy[lev]).array(mfi), (*a_fz[lev]).array(mfi),
                                                    tforces.array(mfi), divu.array(mfi),
                                                    state_in[lev]->array(mfi),state_comp,ncomp,
                                                    u_mac[lev]->array(mfi), v_mac[lev]->array(mfi), w_mac[lev]->array(mfi), bc);

                    }
                    else
                       incflo_compute_fluxes_on_box(lev, bx, (*a_fx[lev])[mfi], (*a_fy[lev])[mfi], (*a_fz[lev])[mfi],
                                                    (*state_in[lev])[mfi], state_comp, ncomp,
                                                    (*xslopes_in[lev])[mfi], (*yslopes_in[lev])[mfi], (*zslopes_in[lev])[mfi], slopes_comp,
                                                    (*u_mac[lev])[mfi], (*v_mac[lev])[mfi], (*w_mac[lev])[mfi]);
#endif
        }// MFIter
}

void
incflo::incflo_compute_fluxes_on_box(const int lev, Box& bx,
                                     FArrayBox& a_fx,
                                     FArrayBox& a_fy,
                                     FArrayBox& a_fz,
                                     const FArrayBox& state_in,
                                     const int state_comp, const int ncomp,
                                     const FArrayBox& xslopes_in,
                                     const FArrayBox& yslopes_in,
                                     const FArrayBox& zslopes_in,
                                     const int slopes_comp,
                                     const FArrayBox& u_mac,
                                     const FArrayBox& v_mac,
                                     const FArrayBox& w_mac)
{
  Box domain(geom[lev].Domain());

  const amrex::Dim3 dom_low = amrex::lbound(domain);
  const amrex::Dim3 dom_high = amrex::ubound(domain);

  Array4<Real> const& fx = a_fx.array();
  Array4<Real> const& fy = a_fy.array();
  Array4<Real> const& fz = a_fz.array();

  Array4<const Real> const& state = state_in.array();

  Array4<const Real> const& u = u_mac.array();
  Array4<const Real> const& v = v_mac.array();
  Array4<const Real> const& w = w_mac.array();

  Array4<const Real> const& x_slopes = xslopes_in.array();
  Array4<const Real> const& y_slopes = yslopes_in.array();
  Array4<const Real> const& z_slopes = zslopes_in.array();

  Array4<int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<int> const& bct_klo = bc_klo[lev]->array();
  Array4<int> const& bct_khi = bc_khi[lev]->array();

  const Box ubx       = amrex::surroundingNodes(bx,0);
  const Box vbx       = amrex::surroundingNodes(bx,1);
  const Box wbx       = amrex::surroundingNodes(bx,2);

  // Vectorize the boundary conditions list in order to use it in lambda
  // functions
  const GpuArray<int, 3> bc_types =
    {bc_list.get_minf(), bc_list.get_pinf(), bc_list.get_pout()};

  amrex::ParallelFor(ubx,ncomp,
    [slopes_comp,state_comp,dom_low,dom_high,bct_ilo,bct_ihi,bc_types,state,x_slopes,u,fx] 
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    Real state_w(0)  ;
    Real state_mns(0); Real state_pls(0);

    //
    // West face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if ((i == dom_low.x) and
     ugradu_aux::is_equal_to_any(bct_ilo(dom_low.x-1,j,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_w = state(i-1,j,k,state_comp+n);

    } else if ((i == dom_high.x+1) and
     ugradu_aux::is_equal_to_any(bct_ihi(dom_high.x+1,j,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_w = state(i,j,k,state_comp+n);
    } else {
      state_pls = state(i  ,j,k,state_comp+n) - .5*x_slopes(i  ,j,k,slopes_comp+n);
      state_mns = state(i-1,j,k,state_comp+n) + .5*x_slopes(i-1,j,k,slopes_comp+n);
      state_w = upwind( state_mns, state_pls, u(i,j,k) );
    }
    fx(i,j,k,n) = u(i,j,k) * state_w;
  });

  amrex::ParallelFor(vbx,ncomp,
    [slopes_comp,state_comp,dom_low,dom_high,bct_jlo,bct_jhi,bc_types,state,y_slopes,v,fy] 
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    Real state_s(0)  ;
    Real state_mns(0); Real state_pls(0);

    //
    // South face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((j == dom_low.y) and
     ugradu_aux::is_equal_to_any(bct_jlo(i,dom_low.y-1,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_s = state(i,j-1,k,state_comp+n);
    } else if ((j == dom_high.y+1) and
     ugradu_aux::is_equal_to_any(bct_jhi(i,dom_high.y+1,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_s = state(i,j,k,state_comp+n);
    } else {
      state_pls = state(i,j  ,k,state_comp+n) - .5*y_slopes(i,j  ,k,slopes_comp+n);
      state_mns = state(i,j-1,k,state_comp+n) + .5*y_slopes(i,j-1,k,slopes_comp+n);

      state_s = upwind( state_mns, state_pls, v(i,j,k) );
    }
    fy(i,j,k,n) = v(i,j,k) * state_s;
  });

  amrex::ParallelFor(wbx,ncomp,
    [slopes_comp,state_comp,dom_low,dom_high,bct_klo,bct_khi,bc_types,state,z_slopes,w,fz] 
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    Real state_b(0)  ;
    Real state_mns(0); Real state_pls(0);
    //
    // Bottom face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((k == dom_low.z) and
     ugradu_aux::is_equal_to_any(bct_klo(i,j,dom_low.z-1,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_b = state(i,j,k-1,state_comp+n);
    } else if ((k == dom_high.z+1) and
     ugradu_aux::is_equal_to_any(bct_khi(i,j,dom_high.z+1,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_b = state(i,j,k,state_comp+n);
    } else {
      state_pls = state(i,j,k  ,state_comp+n) - .5*z_slopes(i,j,k  ,slopes_comp+n);
      state_mns = state(i,j,k-1,state_comp+n) + .5*z_slopes(i,j,k-1,slopes_comp+n);

      state_b = upwind( state_mns, state_pls, w(i,j,k) );
    }
    fz(i,j,k,n) = w(i,j,k) * state_b;
  });
}


#ifdef AMREX_USE_EB
//
// Compute the three components of the convection term when we have embedded
// boundaries
//
void
incflo::incflo_compute_eb_fluxes_on_box(const int lev, Box& bx,
                                        FArrayBox& a_fx,
                                        FArrayBox& a_fy,
                                        FArrayBox& a_fz,
                                        const FArrayBox& state_in,
                                        const int state_comp, const int ncomp,
                                        const FArrayBox& xslopes_in,
                                        const FArrayBox& yslopes_in,
                                        const FArrayBox& zslopes_in,
                                        const int slopes_comp,
                                        const FArrayBox& u_mac,
                                        const FArrayBox& v_mac,
                                        const FArrayBox& w_mac,
                                        const FArrayBox& afrac_x_fab,
                                        const FArrayBox& afrac_y_fab,
                                        const FArrayBox& afrac_z_fab,
                                        const FArrayBox& face_centroid_x,
                                        const FArrayBox& face_centroid_y,
                                        const FArrayBox& face_centroid_z,
                                        const FArrayBox& volfrac,
                                        const FArrayBox& bndry_centroid,
                                        const IArrayBox& cc_mask,
                                        const EBCellFlagFab& flags)
{
  Box domain(geom[lev].Domain());

  const amrex::Dim3 dom_low = amrex::lbound(domain);
  const amrex::Dim3 dom_high = amrex::ubound(domain);

  Array4<Real> const& fx = a_fx.array();
  Array4<Real> const& fy = a_fy.array();
  Array4<Real> const& fz = a_fz.array();

  Array4<const Real> const& state = state_in.array();

  Array4<const Real> const& areafrac_x = afrac_x_fab.array();
  Array4<const Real> const& areafrac_y = afrac_y_fab.array();
  Array4<const Real> const& areafrac_z = afrac_z_fab.array();

  Array4<const Real> const& u = u_mac.array();
  Array4<const Real> const& v = v_mac.array();
  Array4<const Real> const& w = w_mac.array();

  Array4<const Real> const& x_slopes = xslopes_in.array();
  Array4<const Real> const& y_slopes = yslopes_in.array();
  Array4<const Real> const& z_slopes = zslopes_in.array();

  Array4<int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<int> const& bct_klo = bc_klo[lev]->array();
  Array4<int> const& bct_khi = bc_khi[lev]->array();

  const Box ubx       = amrex::surroundingNodes(bx,0);
  const Box vbx       = amrex::surroundingNodes(bx,1);
  const Box wbx       = amrex::surroundingNodes(bx,2);

  const Box ubx_grown = amrex::surroundingNodes(amrex::grow(bx,1),0);
  const Box vbx_grown = amrex::surroundingNodes(amrex::grow(bx,1),1);
  const Box wbx_grown = amrex::surroundingNodes(amrex::grow(bx,1),2);

  FArrayBox s_on_x_face(ubx_grown, ncomp);
  FArrayBox s_on_y_face(vbx_grown, ncomp);
  FArrayBox s_on_z_face(wbx_grown, ncomp);

  // These lines ensure that the temporary Fabs above aren't destroyed
  //   before we're done with them when running with GPUs
  Elixir eli_x = s_on_x_face.elixir();
  Elixir eli_y = s_on_y_face.elixir();
  Elixir eli_z = s_on_z_face.elixir();

  Array4<Real> const& sx = s_on_x_face.array();
  Array4<Real> const& sy = s_on_y_face.array();
  Array4<Real> const& sz = s_on_z_face.array();

  // Face centroids
  const auto& fcx_fab = face_centroid_x.array();
  const auto& fcy_fab = face_centroid_y.array();
  const auto& fcz_fab = face_centroid_z.array();

  const auto& ccm_fab = cc_mask.const_array();

  const GpuArray<int, 3> bc_types =
    {bc_list.get_minf(), bc_list.get_pinf(), bc_list.get_pout()};

  const Real my_huge = 1.2345e300;
  //
  // First compute the convective fluxes at the face center
  // Do this on ALL faces on the tile, i.e. INCLUDE as many ghost faces as
  // possible
  //

  //
  // ===================== X =====================
  //
  amrex::ParallelFor(ubx_grown,ncomp,
    [my_huge,slopes_comp,state_comp,dom_low,dom_high,bct_ilo,bct_ihi,bc_types,areafrac_x,x_slopes,state,u,sx] 
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    Real upls(0); Real umns(0);

    if( areafrac_x(i,j,k) > 0 ) {
      if( i <= dom_low.x and
       ugradu_aux::is_equal_to_any(bct_ilo(dom_low.x-1,j,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        sx(i,j,k,n) = state(dom_low.x-1,j,k,state_comp+n);
      }
      else if( i >= dom_high.x+1 and
       ugradu_aux::is_equal_to_any(bct_ihi(dom_high.x+1,j,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        sx(i,j,k,n) = state(dom_high.x+1,j,k,state_comp+n);
      }
      else {
        upls = state(i  ,j,k,state_comp+n) - .5*x_slopes(i  ,j,k,slopes_comp+n);
        umns = state(i-1,j,k,state_comp+n) + .5*x_slopes(i-1,j,k,slopes_comp+n);

        sx(i,j,k,n) = upwind( umns, upls, u(i,j,k) );
      }
    } else {
        sx(i,j,k,n) = my_huge;
    }
  });

  amrex::ParallelFor(ubx,ncomp, [my_huge,fcx_fab,ccm_fab,areafrac_x,sx,u,fx] 
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    if( areafrac_x(i,j,k) > 0 ) {
       int jj = j + static_cast<int>(std::copysign(1.0, fcx_fab(i,j,k,0)));
       int kk = k + static_cast<int>(std::copysign(1.0, fcx_fab(i,j,k,1)));

       Real fracy = (ccm_fab(i-1,jj,k) || ccm_fab(i,jj,k)) ? std::abs(fcx_fab(i,j,k,0)) : 0.0;
       Real fracz = (ccm_fab(i-1,j,kk) || ccm_fab(i,j,kk)) ? std::abs(fcx_fab(i,j,k,1)) : 0.0;

       Real s_on_x_centroid = (1.0-fracy)*(1.0-fracz)*sx(i, j,k ,n)+
                                   fracy *(1.0-fracz)*sx(i,jj,k ,n)+
                                   fracz *(1.0-fracy)*sx(i, j,kk,n)+
                                   fracy *     fracz *sx(i,jj,kk,n);

       fx(i,j,k,n) = u(i,j,k) * s_on_x_centroid;
    } else {
       fx(i,j,k,n) = my_huge;
    }
  });

  //
  // ===================== Y =====================
  //
  amrex::ParallelFor(vbx_grown,ncomp,
    [my_huge,slopes_comp,state_comp,dom_low,dom_high,bct_jlo,bct_jhi,bc_types,areafrac_y,y_slopes,state,v,sy] 
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    Real vpls(0); Real vmns(0);

    if( areafrac_y(i,j,k) > 0 ) {
      if( j <= dom_low.y and
       ugradu_aux::is_equal_to_any(bct_jlo(i,dom_low.y-1,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        sy(i,j,k,n) = state(i,dom_low.y-1,k,state_comp+n);
      }
      else if( j >= dom_high.y+1 and
       ugradu_aux::is_equal_to_any(bct_jhi(i,dom_high.y+1,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        sy(i,j,k,n) = state(i,dom_high.y+1,k,state_comp+n);
      }
      else {
        vpls = state(i,j  ,k,state_comp+n) - .5*y_slopes(i,j  ,k,slopes_comp+n);
        vmns = state(i,j-1,k,state_comp+n) + .5*y_slopes(i,j-1,k,slopes_comp+n);

        sy(i,j,k,n) = upwind( vmns, vpls, v(i,j,k) );
      }
    }
    else {
        sy(i,j,k,n) = my_huge;
    }
  });

  amrex::ParallelFor(vbx,ncomp, [my_huge,fcy_fab,ccm_fab,areafrac_y,sy,v,fy] 
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    if ( areafrac_y(i,j,k) > 0 ) {
       int ii = i + static_cast<int>(std::copysign(1.0,fcy_fab(i,j,k,0)));
       int kk = k + static_cast<int>(std::copysign(1.0,fcy_fab(i,j,k,1)));

       Real fracx = (ccm_fab(ii,j-1,k) || ccm_fab(ii,j,k)) ? std::abs(fcy_fab(i,j,k,0)) : 0.0;
       Real fracz = (ccm_fab(i,j-1,kk) || ccm_fab(i,j,kk)) ? std::abs(fcy_fab(i,j,k,1)) : 0.0;

       Real s_on_y_centroid = (1.0-fracx)*(1.0-fracz)*sy(i ,j,k ,n)+
                                   fracx *(1.0-fracz)*sy(ii,j,k ,n)+
                                   fracz *(1.0-fracx)*sy(i ,j,kk,n)+
                                   fracx *     fracz *sy(ii,j,kk,n);
       fy(i,j,k,n) = v(i,j,k) * s_on_y_centroid;
    } else
       fy(i,j,k,n) = my_huge;

  });

  //
  // ===================== Z =====================
  //
  amrex::ParallelFor(wbx_grown,ncomp,
    [my_huge,slopes_comp,state_comp,dom_low,dom_high,bct_klo,bct_khi,bc_types,areafrac_z,z_slopes,state,w,sz] 
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    Real wpls(0); Real wmns(0);

    if( areafrac_z(i,j,k) > 0 ) {
      if( k <= dom_low.z and
       ugradu_aux::is_equal_to_any(bct_klo(i,j,dom_low.z-1,0),
                                   bc_types.data(), bc_types.size()))
      {
        sz(i,j,k,n) = state(i,j,dom_low.z-1,state_comp+n);
      }
      else if( k >= dom_high.z+1 and
       ugradu_aux::is_equal_to_any(bct_khi(i,j,dom_high.z+1,0),
                                   bc_types.data(), bc_types.size()))
      {
        sz(i,j,k,n) = state(i,j,dom_high.z+1,state_comp+n);
      }
      else {
        wpls = state(i,j,k  ,state_comp+n) - .5*z_slopes(i,j,k  ,slopes_comp+n);
        wmns = state(i,j,k-1,state_comp+n) + .5*z_slopes(i,j,k-1,slopes_comp+n);

        sz(i,j,k,n) = upwind( wmns, wpls, w(i,j,k) );
      }
    }
    else {
        sz(i,j,k,n) = my_huge;
    }
  });

  amrex::ParallelFor(wbx,ncomp, [my_huge,fcz_fab,ccm_fab,areafrac_z,sz,w,fz] 
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    if( areafrac_z(i,j,k) > 0 ) {
       int ii = i + static_cast<int>(std::copysign(1.0,fcz_fab(i,j,k,0)));
       int jj = j + static_cast<int>(std::copysign(1.0,fcz_fab(i,j,k,1)));

       Real fracx = (ccm_fab(ii,j,k-1) || ccm_fab(ii,j,k)) ? std::abs(fcz_fab(i,j,k,0)) : 0.0;
       Real fracy = (ccm_fab(i,jj,k-1) || ccm_fab(i,jj,k)) ? std::abs(fcz_fab(i,j,k,1)) : 0.0;

       Real s_on_z_centroid = (1.0-fracx)*(1.0-fracy)*sz(i ,j ,k,n)+
                                   fracx *(1.0-fracy)*sz(ii,j ,k,n)+
                                   fracy *(1.0-fracx)*sz(i ,jj,k,n)+
                                   fracx *     fracy *sz(ii,jj,k,n);

       fz(i,j,k,n) = w(i,j,k) * s_on_z_centroid;
    } else
       fz(i,j,k,n) = my_huge;
  });
}
#endif
