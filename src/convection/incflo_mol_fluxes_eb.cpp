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
                                           Array4<Real const> const& ccc)
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

    // ****************************************************************************
    // Predict to x-faces
    // ****************************************************************************

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    auto extdir_lohi = has_extdir(h_bcrec, ncomp, static_cast<int>(Direction::x));
    bool has_extdir_lo = extdir_lohi.first;
    bool has_extdir_hi = extdir_lohi.second;

    if ((has_extdir_lo and domain_ilo >= xbx.smallEnd(0)-1) or
        (has_extdir_hi and domain_ihi <= xbx.bigEnd(0)))
    {
        amrex::ParallelFor(xbx, ncomp,
        [d_bcrec,domain_ilo,domain_jlo,domain_klo,domain_ihi,domain_jhi,domain_khi,
        q,ccc,fcx,flag,umac,small,fx]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           bool extdir_ilo = d_bcrec[n].lo(0) == BCType::ext_dir;
           bool extdir_ihi = d_bcrec[n].hi(0) == BCType::ext_dir;
           bool extdir_jlo = d_bcrec[n].lo(1) == BCType::ext_dir;
           bool extdir_jhi = d_bcrec[n].hi(1) == BCType::ext_dir;
           bool extdir_klo = d_bcrec[n].lo(2) == BCType::ext_dir;
           bool extdir_khi = d_bcrec[n].hi(2) == BCType::ext_dir;
           Real qs;

           if (flag(i,j,k).isConnected(-1,0,0)) 
           {
               if (extdir_ilo and i <= domain_ilo) {
                   qs = q(domain_ilo-1,j,k,n);
               } else if (extdir_ihi and i >= domain_ihi+1) {
                   qs = q(domain_ihi+1,j,k,n);
               } else {

                   Real yf = fcx(i,j,k,0); // local (y,z) of centroid of z-face we are extrapolating to
                   Real zf = fcx(i,j,k,1);
 
                   Real xc = ccc(i,j,k,0); // centroid of cell (i,j,k)
                   Real yc = ccc(i,j,k,1);
                   Real zc = ccc(i,j,k,2);
 
                   Real delta_x = 0.5 + xc;
                   Real delta_y = yf  - yc;
                   Real delta_z = zf  - zc;
    
                   Real cc_qmax = std::max(q(i,j,k,n),q(i-1,j,k,n));
                   Real cc_qmin = std::min(q(i,j,k,n),q(i-1,j,k,n));
     
                   // Compute slopes of component "n" of q
                   const auto& slopes_eb_hi = incflo_slopes_extdir_eb(i,j,k,n,q,ccc,flag,
                                              extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                              extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                              extdir_klo, extdir_khi, domain_klo, domain_khi);
 
                   Real qpls = q(i  ,j,k,n) - delta_x * slopes_eb_hi[0]
                                            + delta_y * slopes_eb_hi[1]
                                            + delta_z * slopes_eb_hi[2];
    
                   qpls = std::max(std::min(qpls, cc_qmax), cc_qmin);
    
                   xc = ccc(i-1,j,k,0); // centroid of cell (i-1,j,k)
                   yc = ccc(i-1,j,k,1);
                   zc = ccc(i-1,j,k,2);
    
                   delta_x = 0.5 - xc;
                   delta_y = yf  - yc;
                   delta_z = zf  - zc;

                   // Compute slopes of component "n" of q
                   const auto& slopes_eb_lo = incflo_slopes_extdir_eb(i-1,j,k,n,q,ccc,flag,
                                              extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                              extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                              extdir_klo, extdir_khi, domain_klo, domain_khi);
    
                   Real qmns = q(i-1,j,k,n) + delta_x * slopes_eb_lo[0]
                                            + delta_y * slopes_eb_lo[1]
                                            + delta_z * slopes_eb_lo[2];
    
                   qmns = std::max(std::min(qmns, cc_qmax), cc_qmin);

                    if (umac(i,j,k) > small) {
                        qs = qmns;
                    } else if (umac(i,j,k) < -small) {
                        qs = qpls;
                    } else {
                        qs = 0.5*(qmns+qpls);
                    }
               }

               fx(i,j,k,n) = umac(i,j,k) * qs;
   
           } else {
               fx(i,j,k,n) = 0.0;
           }
        });
    }
    else
    {
        amrex::ParallelFor(xbx, ncomp,
        [q,ccc,fcx,flag,umac,small,fx]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           Real qs;
           if (flag(i,j,k).isConnected(-1,0,0)) 
           {
               Real yf = fcx(i,j,k,0); // local (y,z) of centroid of z-face we are extrapolating to
               Real zf = fcx(i,j,k,1);
 
               Real xc = ccc(i,j,k,0); // centroid of cell (i,j,k)
               Real yc = ccc(i,j,k,1);
               Real zc = ccc(i,j,k,2);

               Real delta_x = 0.5 + xc;
               Real delta_y = yf  - yc;
               Real delta_z = zf  - zc;

               Real cc_qmax = std::max(q(i,j,k,n),q(i-1,j,k,n));
               Real cc_qmin = std::min(q(i,j,k,n),q(i-1,j,k,n));

               // Compute slopes of component "n" of q
               const auto& slopes_eb_hi = incflo_slopes_eb(i,j,k,n,q,ccc,flag);

               Real qpls = q(i  ,j,k,n) - delta_x * slopes_eb_hi[0]
                                        + delta_y * slopes_eb_hi[1]
                                        + delta_z * slopes_eb_hi[2];

               qpls = std::max(std::min(qpls, cc_qmax), cc_qmin);

               xc = ccc(i-1,j,k,0); // centroid of cell (i-1,j,k)
               yc = ccc(i-1,j,k,1);
               zc = ccc(i-1,j,k,2);

               delta_x = 0.5 - xc;
               delta_y = yf  - yc;
               delta_z = zf  - zc;

               // Compute slopes of component "n" of q
               const auto& slopes_eb_lo = incflo_slopes_eb(i-1,j,k,n,q,ccc,flag);

               Real qmns = q(i-1,j,k,n) + delta_x * slopes_eb_lo[0]
                                        + delta_y * slopes_eb_lo[1]
                                        + delta_z * slopes_eb_lo[2];
    
               qmns = std::max(std::min(qmns, cc_qmax), cc_qmin);

               if (umac(i,j,k) > small) {
                    qs = qmns;
                } else if (umac(i,j,k) < -small) {
                    qs = qpls;
                } else {
                    qs = 0.5*(qmns+qpls);
                }

                fx(i,j,k,n) = umac(i,j,k) * qs;

           } else {
               fx(i,j,k,n) = 0.0;
           }
        });
    }

    // ****************************************************************************
    // Predict to y-faces
    // ****************************************************************************

    extdir_lohi = has_extdir(h_bcrec, ncomp,  static_cast<int>(Direction::y));
    has_extdir_lo = extdir_lohi.first;
    has_extdir_hi = extdir_lohi.second;
    if ((has_extdir_lo and domain_jlo >= ybx.smallEnd(1)-1) or
        (has_extdir_hi and domain_jhi <= ybx.bigEnd(1)))
    {
        amrex::ParallelFor(ybx, ncomp,
        [d_bcrec,domain_ilo,domain_jlo,domain_klo,domain_ihi,domain_jhi,domain_khi,
         q,ccc,fcy,flag,vmac,small,fy]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qs;
            if (flag(i,j,k).isConnected(0,-1,0)) 
            {
                bool extdir_ilo = d_bcrec[n].lo(0) == BCType::ext_dir;
                bool extdir_ihi = d_bcrec[n].hi(0) == BCType::ext_dir;
                bool extdir_jlo = d_bcrec[n].lo(1) == BCType::ext_dir;
                bool extdir_jhi = d_bcrec[n].hi(1) == BCType::ext_dir;
                bool extdir_klo = d_bcrec[n].lo(2) == BCType::ext_dir;
                bool extdir_khi = d_bcrec[n].hi(2) == BCType::ext_dir;

                if (extdir_jlo and j <= domain_jlo) {
                    qs = q(i,domain_jlo-1,k,n);
                } else if (extdir_jhi and j >= domain_jhi+1) {
                    qs = q(i,domain_jhi+1,k,n);
                } else {

                   Real xf = fcy(i,j,k,0); // local (x,z) of centroid of z-face we are extrapolating to
                   Real zf = fcy(i,j,k,1);
 
                   Real xc = ccc(i,j,k,0); // centroid of cell (i,j,k)
                   Real yc = ccc(i,j,k,1);
                   Real zc = ccc(i,j,k,2);
 
                   Real delta_x = xf  - xc;
                   Real delta_y = 0.5 + yc;
                   Real delta_z = zf  - zc;
    
                   Real cc_qmax = std::max(q(i,j,k,n),q(i,j-1,k,n));
                   Real cc_qmin = std::min(q(i,j,k,n),q(i,j-1,k,n));
     
                   // Compute slopes of component "n" of q
                   const auto& slopes_eb_hi = incflo_slopes_extdir_eb(i,j,k,n,q,ccc,flag,
                                              extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                              extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                              extdir_klo, extdir_khi, domain_klo, domain_khi);
 
                   Real qpls = q(i,j  ,k,n) + delta_x * slopes_eb_hi[0]
                                            - delta_y * slopes_eb_hi[1]
                                            + delta_z * slopes_eb_hi[2];
    
                   qpls = std::max(std::min(qpls, cc_qmax), cc_qmin);
    
                   xc = ccc(i,j-1,k,0); // centroid of cell (i-1,j,k)
                   yc = ccc(i,j-1,k,1);
                   zc = ccc(i,j-1,k,2);
    
                   delta_x = xf  - xc;
                   delta_y = 0.5 - yc;
                   delta_z = zf  - zc;

                   // Compute slopes of component "n" of q
                   const auto& slopes_eb_lo = incflo_slopes_extdir_eb(i,j-1,k,n,q,ccc,flag,
                                              extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                              extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                              extdir_klo, extdir_khi, domain_klo, domain_khi);
    
                   Real qmns = q(i,j-1,k,n) + delta_x * slopes_eb_lo[0]
                                            + delta_y * slopes_eb_lo[1]
                                            + delta_z * slopes_eb_lo[2];
    
                   qmns = std::max(std::min(qmns, cc_qmax), cc_qmin);

                    if (vmac(i,j,k) > small) {
                        qs = qmns;
                    } else if (vmac(i,j,k) < -small) {
                        qs = qpls;
                    } else {
                        qs = 0.5*(qmns+qpls);
                    }
                }

                fy(i,j,k,n) = vmac(i,j,k) * qs;

           } else {
                fy(i,j,k,n) = 0.0;
           }
        });
    }
    else
    {
        amrex::ParallelFor(ybx, ncomp,
        [q,ccc,fcy,flag,vmac,small,fy]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real qs;
            if (flag(i,j,k).isConnected(0,-1,0)) {

               Real xf = fcy(i,j,k,0); // local (x,z) of centroid of z-face we are extrapolating to
               Real zf = fcy(i,j,k,1);
 
               Real xc = ccc(i,j,k,0); // centroid of cell (i,j,k)
               Real yc = ccc(i,j,k,1);
               Real zc = ccc(i,j,k,2);

               Real delta_x = xf  - xc;
               Real delta_y = 0.5 + yc;
               Real delta_z = zf  - zc;

               Real cc_qmax = std::max(q(i,j,k,n),q(i,j-1,k,n));
               Real cc_qmin = std::min(q(i,j,k,n),q(i,j-1,k,n));

               // Compute slopes of component "n" of q
               const auto& slopes_eb_hi = incflo_slopes_eb(i,j,k,n,q,ccc,flag);

               Real qpls = q(i,j  ,k,n) + delta_x * slopes_eb_hi[0]
                                        - delta_y * slopes_eb_hi[1]
                                        + delta_z * slopes_eb_hi[2];

               qpls = std::max(std::min(qpls, cc_qmax), cc_qmin);

               xc = ccc(i,j-1,k,0); // centroid of cell (i-1,j,k)
               yc = ccc(i,j-1,k,1);
               zc = ccc(i,j-1,k,2);

               delta_x = xf  - xc;
               delta_y = 0.5 - yc;
               delta_z = zf  - zc;

               // Compute slopes of component "n" of q
               const auto& slopes_eb_lo = incflo_slopes_eb(i,j-1,k,n,q,ccc,flag);

               Real qmns = q(i,j-1,k,n) + delta_x * slopes_eb_lo[0]
                                        + delta_y * slopes_eb_lo[1]
                                        + delta_z * slopes_eb_lo[2];
    
               qmns = std::max(std::min(qmns, cc_qmax), cc_qmin);

               if (vmac(i,j,k) > small) {
                   qs = qmns;
               } else if (vmac(i,j,k) < -small) {
                   qs = qpls;
               } else {
                   qs = 0.5*(qmns+qpls);
               }

               fy(i,j,k,n) = vmac(i,j,k) * qs;

           } else {
               fy(i,j,k,n) = 0.0;
           }
        });
    }

    // ****************************************************************************
    // Predict to z-faces
    // ****************************************************************************

    extdir_lohi = has_extdir(h_bcrec, ncomp, static_cast<int>(Direction::z));
    has_extdir_lo = extdir_lohi.first;
    has_extdir_hi = extdir_lohi.second;
    if ((has_extdir_lo and domain_klo >= zbx.smallEnd(2)-1) or
        (has_extdir_hi and domain_khi <= zbx.bigEnd(2)))
    {
        amrex::ParallelFor(zbx, ncomp,
        [d_bcrec,domain_ilo,domain_jlo,domain_klo,domain_ihi,domain_jhi,domain_khi,
         q,ccc,fcz,flag,wmac,small,fz]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (flag(i,j,k).isConnected(0,0,-1)) {

                bool extdir_ilo = d_bcrec[n].lo(0) == BCType::ext_dir;
                bool extdir_ihi = d_bcrec[n].hi(0) == BCType::ext_dir;
                bool extdir_jlo = d_bcrec[n].lo(1) == BCType::ext_dir;
                bool extdir_jhi = d_bcrec[n].hi(1) == BCType::ext_dir;
                bool extdir_klo = d_bcrec[n].lo(2) == BCType::ext_dir;
                bool extdir_khi = d_bcrec[n].hi(2) == BCType::ext_dir;

                Real qs;
                if (extdir_klo and k <= domain_klo) {
                    qs = q(i,j,domain_klo-1,n);
                } else if (extdir_khi and k >= domain_khi+1) {
                    qs = q(i,j,domain_khi+1,n);
                } else {

                    Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
                    Real yf = fcz(i,j,k,1);
 
                    Real xc = ccc(i,j,k,0); // centroid of cell (i,j,k)
                    Real yc = ccc(i,j,k,1);
                    Real zc = ccc(i,j,k,2);
 
                    Real delta_x = xf  - xc;
                    Real delta_y = yf  - yc;
                    Real delta_z = 0.5 + zc;
     
                    Real cc_qmax = std::max(q(i,j,k,n),q(i,j,k-1,n));
                    Real cc_qmin = std::min(q(i,j,k,n),q(i,j,k-1,n));
     
                    // Compute slopes of component "n" of q
                    const auto& slopes_eb_hi = incflo_slopes_extdir_eb(i,j,k,n,q,ccc,flag,
                                               extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                               extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                               extdir_klo, extdir_khi, domain_klo, domain_khi);
 
                    Real qpls = q(i,j,k  ,n) + delta_x * slopes_eb_hi[0]
                                             + delta_y * slopes_eb_hi[1]
                                             - delta_z * slopes_eb_hi[2];
     
                    qpls = std::max(std::min(qpls, cc_qmax), cc_qmin);
     
                    xc = ccc(i,j,k-1,0); // centroid of cell (i,j,k-1)
                    yc = ccc(i,j,k-1,1);
                    zc = ccc(i,j,k-1,2);
     
                    delta_x = xf  - xc;
                    delta_y = yf  - yc;
                    delta_z = 0.5 - zc;

                    // Compute slopes of component "n" of q
                    const auto& slopes_eb_lo = incflo_slopes_extdir_eb(i,j,k-1,n,q,ccc,flag,
                                               extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                               extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                               extdir_klo, extdir_khi, domain_klo, domain_khi);
    
                    Real qmns = q(i,j,k-1,n) + delta_x * slopes_eb_lo[0]
                                             + delta_y * slopes_eb_lo[1]
                                             + delta_z * slopes_eb_lo[2];
    
                    qmns = std::max(std::min(qmns, cc_qmax), cc_qmin);

                    if (wmac(i,j,k) > small) {
                        qs = qmns;
                    } else if (wmac(i,j,k) < -small) {
                        qs = qpls;
                    } else {
                        qs = 0.5*(qmns+qpls);
                    }
                }

                fz(i,j,k,n) = wmac(i,j,k) * qs;

           } else {
                fz(i,j,k,n) = 0.0;
           }
        });
    }
    else
    {
        amrex::ParallelFor(zbx, ncomp,
        [q,ccc,fcz,flag,wmac,small,fz]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (flag(i,j,k).isConnected(0,0,-1)) {
                Real qs;

                Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
                Real yf = fcz(i,j,k,1);
 
                Real xc = ccc(i,j,k,0); // centroid of cell (i,j,k)
                Real yc = ccc(i,j,k,1);
                Real zc = ccc(i,j,k,2);
 
                Real delta_x = xf  - xc;
                Real delta_y = yf  - yc;
                Real delta_z = 0.5 + zc;
     
                Real cc_qmax = std::max(q(i,j,k,n),q(i,j,k-1,n));
                Real cc_qmin = std::min(q(i,j,k,n),q(i,j,k-1,n));
     
                // Compute slopes of component "n" of q
                const auto& slopes_eb_hi = incflo_slopes_eb(i,j,k,n,q,ccc,flag);
 
                Real qpls = q(i,j,k  ,n) + delta_x * slopes_eb_hi[0]
                                         + delta_y * slopes_eb_hi[1]
                                         - delta_z * slopes_eb_hi[2];
 
                qpls = std::max(std::min(qpls, cc_qmax), cc_qmin);
 
                xc = ccc(i,j,k-1,0); // centroid of cell (i,j,k-1)
                yc = ccc(i,j,k-1,1);
                zc = ccc(i,j,k-1,2);
 
                delta_x = xf  - xc;
                delta_y = yf  - yc;
                delta_z = 0.5 - zc;

                // Compute slopes of component "n" of q
                const auto& slopes_eb_lo = incflo_slopes_eb(i,j,k-1,n,q,ccc,flag);
    
                Real qmns = q(i,j,k-1,n) + delta_x * slopes_eb_lo[0]
                                         + delta_y * slopes_eb_lo[1]
                                         + delta_z * slopes_eb_lo[2];

                qmns = std::max(std::min(qmns, cc_qmax), cc_qmin);

                if (wmac(i,j,k) > small) {
                    qs = qmns;
                } else if (wmac(i,j,k) < -small) {
                    qs = qpls;
                } else {
                    qs = 0.5*(qmns+qpls);
                }

                fz(i,j,k,n) = wmac(i,j,k) * qs;

           } else {
                fz(i,j,k,n) = 0.0;
           }
        });
    }
}
#endif
