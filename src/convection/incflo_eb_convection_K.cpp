#include <incflo_convection_K.H>
#include <incflo.H>

using namespace amrex;

#ifdef AMREX_USE_EB
void incflo::predict_vels_on_faces_eb (int lev, Box const& ccbx,
                                       Box const& ubx, Box const& vbx, Box const& wbx,
                                       Array4<Real> const& u, Array4<Real> const& v,
                                       Array4<Real> const& w, Array4<Real const> const& vcc,
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

    auto const bc_ilo = m_bc_type[Orientation(Direction::x,Orientation::low)];
    auto const bc_ihi = m_bc_type[Orientation(Direction::x,Orientation::high)];
    auto const bc_jlo = m_bc_type[Orientation(Direction::y,Orientation::low)];
    auto const bc_jhi = m_bc_type[Orientation(Direction::y,Orientation::high)];
    auto const bc_klo = m_bc_type[Orientation(Direction::z,Orientation::low)];
    auto const bc_khi = m_bc_type[Orientation(Direction::z,Orientation::high)];

    bool extdir_ilo = (bc_ilo == BC::mass_inflow) or (bc_ilo == BC::no_slip_wall);
    bool extdir_ihi = (bc_ihi == BC::mass_inflow) or (bc_ihi == BC::no_slip_wall);
    bool extdir_jlo = (bc_jlo == BC::mass_inflow) or (bc_jlo == BC::no_slip_wall);
    bool extdir_jhi = (bc_jhi == BC::mass_inflow) or (bc_jhi == BC::no_slip_wall);
    bool extdir_klo = (bc_klo == BC::mass_inflow) or (bc_klo == BC::no_slip_wall);
    bool extdir_khi = (bc_khi == BC::mass_inflow) or (bc_khi == BC::no_slip_wall);

    // At an ext_dir boundary, the boundary value is on the face, not cell center.

    Real small_vel = 1.e-10;

    // ****************************************************************************
    // Predict to x-faces
    // ****************************************************************************
    if ((extdir_ilo and domain_ilo >= ubx.smallEnd(0)-1) or
        (extdir_ihi and domain_ihi <= ubx.bigEnd(0)))
    {
        amrex::ParallelFor(Box(ubx),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (flag(i,j,k).isConnected(-1,0,0))
            {
               Real yf = fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
               Real zf = fcx(i,j,k,1);

               Real xc = ccc(i,j,k,0); // centroid of cell (i,j,k)
               Real yc = ccc(i,j,k,1);
               Real zc = ccc(i,j,k,2);

               Real delta_x = 0.5 + xc;
               Real delta_y = yf  - yc;
               Real delta_z = zf  - zc;

               Real cc_umax = std::max(vcc(i,j,k,0), vcc(i-1,j,k,0));
               Real cc_umin = std::min(vcc(i,j,k,0), vcc(i-1,j,k,0));

               // Compute slopes of component "0" of vcc
               const auto& slopes_eb_hi = incflo_slopes_extdir_eb(i,j,k,0,vcc,ccc,flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real upls = vcc(i  ,j,k,0) - delta_x * slopes_eb_hi[0]
                                          + delta_y * slopes_eb_hi[1]
                                          + delta_z * slopes_eb_hi[2];

               upls = std::max(std::min(upls, cc_umax), cc_umin);

               xc = ccc(i-1,j,k,0); // centroid of cell (i-1,j,k)
               yc = ccc(i-1,j,k,1);
               zc = ccc(i-1,j,k,2);

               delta_x = 0.5 - xc;
               delta_y = yf  - yc;
               delta_z = zf  - zc;

               // Compute slopes of component "0" of vcc
               const auto& slopes_eb_lo = incflo_slopes_extdir_eb(i-1,j,k,0,vcc,ccc,flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real umns = vcc(i-1,j,k,0) + delta_x * slopes_eb_lo[0]
                                          + delta_y * slopes_eb_lo[1]
                                          + delta_z * slopes_eb_lo[2];

               umns = std::max(std::min(umns, cc_umax), cc_umin);

               if ( umns < 0.0 && upls > 0.0 ) {
                  u(i,j,k) = 0.0;
               } else {
                  Real avg = 0.5 * ( upls + umns );
                  if ( std::abs(avg) <  small_vel) { u(i,j,k) = 0.0;
                  } else if (avg >= 0)             { u(i,j,k) = umns;
                  } else                           { u(i,j,k) = upls;
                  }
               }

               if (extdir_ilo and i == domain_ilo) {
                   u(i,j,k) = vcc(i-1,j,k,0);
               } else if (extdir_ihi and i == domain_ihi+1) {
                   u(i,j,k) = vcc(i,j,k,0);
               }

            } else {
               u(i,j,k) = 0.0;
            } 
        });
    }
    else
    {
        amrex::ParallelFor(Box(ubx),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (flag(i,j,k).isConnected(-1,0,0))
            {
               Real yf = fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
               Real zf = fcx(i,j,k,1);

               Real xc = ccc(i,j,k,0); // centroid of cell (i,j,k)
               Real yc = ccc(i,j,k,1);
               Real zc = ccc(i,j,k,2);

               Real delta_x = 0.5 + xc;
               Real delta_y = yf  - yc;
               Real delta_z = zf  - zc;

               Real cc_umax = std::max(vcc(i,j,k,0), vcc(i-1,j,k,0));
               Real cc_umin = std::min(vcc(i,j,k,0), vcc(i-1,j,k,0));

               // Compute slopes of component "0" of vcc
               const auto slopes_eb_hi = incflo_slopes_eb(i,j,k,0,vcc,ccc,flag);

               Real upls = vcc(i  ,j,k,0) - delta_x * slopes_eb_hi[0]
                                          + delta_y * slopes_eb_hi[1]
                                          + delta_z * slopes_eb_hi[2];

               upls = std::max(std::min(upls, cc_umax), cc_umin);

               xc = ccc(i-1,j,k,0); // centroid of cell (i-1,j,k)
               yc = ccc(i-1,j,k,1);
               zc = ccc(i-1,j,k,2);

               delta_x = 0.5 - xc;
               delta_y = yf  - yc;
               delta_z = zf  - zc;

               // Compute slopes of component "0" of vcc
               const auto& slopes_eb_lo = incflo_slopes_eb(i-1,j,k,0,vcc,ccc,flag);

               Real umns = vcc(i-1,j,k,0) + delta_x * slopes_eb_lo[0]
                                          + delta_y * slopes_eb_lo[1]
                                          + delta_z * slopes_eb_lo[2];

               umns = std::max(std::min(umns, cc_umax), cc_umin);

               if ( umns < 0.0 && upls > 0.0 ) {
                  u(i,j,k) = 0.0;
               } else {
                  Real avg = 0.5 * ( upls + umns );
                  if ( std::abs(avg) <  small_vel) { u(i,j,k) = 0.0;
                  } else if (avg >= 0)             { u(i,j,k) = umns;
                  } else                           { u(i,j,k) = upls;
                  }
               }

            } else {
               u(i,j,k) = 0.0;
            } 
        });
    }

    // ****************************************************************************
    // Predict to y-faces
    // ****************************************************************************
    if ((extdir_jlo and domain_jlo >= vbx.smallEnd(1)-1) or
        (extdir_jhi and domain_jhi <= vbx.bigEnd(1)))
    {
        amrex::ParallelFor(Box(vbx),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (flag(i,j,k).isConnected(0,-1,0))
            {
               Real xf = fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
               Real zf = fcy(i,j,k,1);

               Real xc = ccc(i,j,k,0); // centroid of cell (i,j,k)
               Real yc = ccc(i,j,k,1);
               Real zc = ccc(i,j,k,2);

               Real delta_x = xf  - xc;
               Real delta_y = 0.5 + yc;
               Real delta_z = zf  - zc;

               Real cc_vmax = std::max(vcc(i,j,k,1), vcc(i,j-1,k,1));
               Real cc_vmin = std::min(vcc(i,j,k,1), vcc(i,j-1,k,1));

               // Compute slopes of component "1" of vcc
               const auto& slopes_eb_hi = incflo_slopes_extdir_eb(i,j,k,1,vcc,ccc,flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real vpls = vcc(i,j  ,k,1) + delta_x * slopes_eb_hi[0]
                                          - delta_y * slopes_eb_hi[1]
                                          + delta_z * slopes_eb_hi[2];

               vpls = std::max(std::min(vpls, cc_vmax), cc_vmin);

               xc = ccc(i,j-1,k,0); // centroid of cell (i,j-1,k)
               yc = ccc(i,j-1,k,1);
               zc = ccc(i,j-1,k,2);

               delta_x = xf  - xc;
               delta_y = 0.5 - yc;
               delta_z = zf  - zc;

               // Compute slopes of component "1" of vcc
               const auto& slopes_eb_lo = incflo_slopes_extdir_eb(i,j-1,k,1,vcc,ccc,flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real vmns = vcc(i,j-1,k,1) + delta_x * slopes_eb_lo[0]
                                          + delta_y * slopes_eb_lo[1]
                                          + delta_z * slopes_eb_lo[2];

               vmns = std::max(std::min(vmns, cc_vmax), cc_vmin);

               if ( vmns < 0.0 && vpls > 0.0 ) {
                  v(i,j,k) = 0.0;
               } else {
                  Real avg = 0.5 * ( vpls + vmns );
                  if ( std::abs(avg) <  small_vel) { v(i,j,k) = 0.0;
                  } else if (avg >= 0)             { v(i,j,k) = vmns;
                  } else                           { v(i,j,k) = vpls;
                  }
               }

               if (extdir_jlo and j == domain_jlo) {
                   v(i,j,k) = vcc(i,j-1,k,1);
               } else if (extdir_jhi and j == domain_jhi+1) {
                   v(i,j,k) = vcc(i,j,k,1);
               }

            } else {
               v(i,j,k) = 0.0;
            } 
        });
    }
    else
    {
        amrex::ParallelFor(Box(vbx),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (flag(i,j,k).isConnected(0,-1,0))
            {
               Real xf = fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
               Real zf = fcy(i,j,k,1);

               Real xc = ccc(i,j,k,0); // centroid of cell (i,j,k)
               Real yc = ccc(i,j,k,1);
               Real zc = ccc(i,j,k,2);

               Real delta_x = xf  - xc;
               Real delta_y = 0.5 + yc;
               Real delta_z = zf  - zc;

               Real cc_vmax = std::max(vcc(i,j,k,1), vcc(i,j-1,k,1));
               Real cc_vmin = std::min(vcc(i,j,k,1), vcc(i,j-1,k,1));

               // Compute slopes of component "1" of vcc
               const auto slopes_eb_hi = incflo_slopes_eb(i,j,k,1,vcc,ccc,flag);

               Real vpls = vcc(i,j  ,k,1) + delta_x * slopes_eb_hi[0]
                                          - delta_y * slopes_eb_hi[1]
                                          + delta_z * slopes_eb_hi[2];

               vpls = std::max(std::min(vpls, cc_vmax), cc_vmin);

               xc = ccc(i,j-1,k,0); // centroid of cell (i,j-1,k)
               yc = ccc(i,j-1,k,1);
               zc = ccc(i,j-1,k,2);

               delta_x = xf  - xc;
               delta_y = 0.5 - yc;
               delta_z = zf  - zc;

               // Compute slopes of component "1" of vcc
               const auto& slopes_eb_lo = incflo_slopes_eb(i,j-1,k,1,vcc,ccc,flag);

               Real vmns = vcc(i,j-1,k,1) + delta_x * slopes_eb_lo[0]
                                          + delta_y * slopes_eb_lo[1]
                                          + delta_z * slopes_eb_lo[2];
                                          
               vmns = std::max(std::min(vmns, cc_vmax), cc_vmin);

               if ( vmns < 0.0 && vpls > 0.0 ) {
                  v(i,j,k) = 0.0;
               } else {
                  Real avg = 0.5 * ( vpls + vmns );
                  if ( std::abs(avg) <  small_vel) { v(i,j,k) = 0.0;
                  } else if (avg >= 0)             { v(i,j,k) = vmns;
                  } else                           { v(i,j,k) = vpls;
                  }
               }

            } else {
               v(i,j,k) = 0.0;
            } 
        });
    }

    // ****************************************************************************
    // Predict to z-faces
    // ****************************************************************************
    if ((extdir_klo and domain_klo >= wbx.smallEnd(2)-1) or
        (extdir_khi and domain_khi <= wbx.bigEnd(2)))
    {
        amrex::ParallelFor(Box(wbx),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (flag(i,j,k).isConnected(0,0,-1))
            {
               Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
               Real yf = fcz(i,j,k,1);

               Real xc = ccc(i,j,k,0); // centroid of cell (i,j,k)
               Real yc = ccc(i,j,k,1);
               Real zc = ccc(i,j,k,2);

               Real delta_x = xf  - xc;
               Real delta_y = yf  - yc;
               Real delta_z = 0.5 + zc;

               Real cc_wmax = std::max(vcc(i,j,k,2), vcc(i,j,k-1,2));
               Real cc_wmin = std::min(vcc(i,j,k,2), vcc(i,j,k-1,2));

               // Compute slopes of component "2" of vcc
               const auto& slopes_eb_hi = incflo_slopes_extdir_eb(i,j,k,2,vcc,ccc,flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real wpls = vcc(i,j,k  ,2) + delta_x * slopes_eb_hi[0]
                                          + delta_y * slopes_eb_hi[1]
                                          - delta_z * slopes_eb_hi[2];

               wpls = std::max(std::min(wpls, cc_wmax), cc_wmin);

               xc = ccc(i,j,k-1,0); // centroid of cell (i,j,k-1)
               yc = ccc(i,j,k-1,1);
               zc = ccc(i,j,k-1,2);

               delta_x = xf  - xc;
               delta_y = yf  - yc;
               delta_z = 0.5 - zc;

               // Compute slopes of component "2" of vcc
               const auto& slopes_eb_lo = incflo_slopes_extdir_eb(i,j,k-1,2,vcc,ccc,flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real wmns = vcc(i,j,k-1,2) + delta_x * slopes_eb_lo[0]
                                          + delta_y * slopes_eb_lo[1]
                                          + delta_z * slopes_eb_lo[2];

               wmns = std::max(std::min(wmns, cc_wmax), cc_wmin);

               if ( wmns < 0.0 && wpls > 0.0 ) {
                  w(i,j,k) = 0.0;
               } else {
                  Real avg = 0.5 * ( wpls + wmns );
                  if ( std::abs(avg) <  small_vel) { w(i,j,k) = 0.0;
                  } else if (avg >= 0)             { w(i,j,k) = wmns;
                  } else                           { w(i,j,k) = wpls;
                  }
               }

                if (extdir_klo and k == domain_klo) {
                    w(i,j,k) = vcc(i,j,k-1,2);
                } else if (extdir_khi and k == domain_khi+1) {
                    w(i,j,k) = vcc(i,j,k,2);
                }

            } else {
               w(i,j,k) = 0.0;
            } 
        });
    }
    else
    {
        amrex::ParallelFor(Box(wbx),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (flag(i,j,k).isConnected(0,0,-1))
            {
               Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
               Real yf = fcz(i,j,k,1);

               Real xc = ccc(i,j,k,0); // centroid of cell (i,j,k)
               Real yc = ccc(i,j,k,1);
               Real zc = ccc(i,j,k,2);

               Real delta_x = xf  - xc;
               Real delta_y = yf  - yc;
               Real delta_z = 0.5 + zc;

               Real cc_wmax = std::max(vcc(i,j,k,2), vcc(i,j,k-1,2));
               Real cc_wmin = std::min(vcc(i,j,k,2), vcc(i,j,k-1,2));

               // Compute slopes of component "2" of vcc
               const auto slopes_eb_hi = incflo_slopes_eb(i,j,k,2,vcc,ccc,flag);

               Real wpls = vcc(i,j,k  ,2) + delta_x * slopes_eb_hi[0]
                                          + delta_y * slopes_eb_hi[1]
                                          - delta_z * slopes_eb_hi[2];

               wpls = std::max(std::min(wpls, cc_wmax), cc_wmin);

               xc = ccc(i,j,k-1,0); // centroid of cell (i,j,k-1)
               yc = ccc(i,j,k-1,1);
               zc = ccc(i,j,k-1,2);

               delta_x = xf  - xc;
               delta_y = yf  - yc;
               delta_z = 0.5 - zc;

               // Compute slopes of component "2" of vcc
               const auto& slopes_eb_lo = incflo_slopes_eb(i,j,k-1,2,vcc,ccc,flag);

               Real wmns = vcc(i,j,k-1,2) + delta_x * slopes_eb_lo[0]
                                          + delta_y * slopes_eb_lo[1]
                                          + delta_z * slopes_eb_lo[2];

               wmns = std::max(std::min(wmns, cc_wmax), cc_wmin);

               if ( wmns < 0.0 && wpls > 0.0 ) {
                  w(i,j,k) = 0.0;
               } else {
                  Real avg = 0.5 * ( wpls + wmns );
                  if ( std::abs(avg) <  small_vel) { w(i,j,k) = 0.0;
                  } else if (avg >= 0)             { w(i,j,k) = wmns;
                  } else                           { w(i,j,k) = wpls;
                  }
               }

            } else {
               w(i,j,k) = 0.0;
            } 
        });
    }
}

void incflo::compute_convective_rate_eb (int lev, Box const& bx, int ncomp,
                                         Array4<Real> const& dUdt,
                                         Array4<Real const> const& fx,
                                         Array4<Real const> const& fy,
                                         Array4<Real const> const& fz,
                                         Array4<EBCellFlag const> const& flag,
                                         Array4<Real const> const& vfrac,
                                         Array4<Real const> const& apx,
                                         Array4<Real const> const& apy,
                                         Array4<Real const> const& apz)
{
    const auto dxinv = Geom(lev).InvCellSizeArray();
    const Box dbox = Geom(lev).growPeriodicDomain(2);
    amrex::ParallelFor(bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (!dbox.contains(IntVect(i,j,k)) or flag(i,j,k).isCovered()) {
            dUdt(i,j,k,n) = 0.0;
        } else if (flag(i,j,k).isRegular()) {
            dUdt(i,j,k,n) = dxinv[0] * (fx(i,j,k,n) - fx(i+1,j,k,n))
                +           dxinv[1] * (fy(i,j,k,n) - fy(i,j+1,k,n))
                +           dxinv[2] * (fz(i,j,k,n) - fz(i,j,k+1,n));
        } else {
            dUdt(i,j,k,n) = (1.0/vfrac(i,j,k)) *
                ( dxinv[0] * (apx(i,j,k)*fx(i,j,k,n) - apx(i+1,j,k)*fx(i+1,j,k,n))
                + dxinv[1] * (apy(i,j,k)*fy(i,j,k,n) - apy(i,j+1,k)*fy(i,j+1,k,n))
                + dxinv[2] * (apz(i,j,k)*fz(i,j,k,n) - apz(i,j,k+1)*fz(i,j,k+1,n)) );
        }
    });
}

void incflo::redistribute_eb (int lev, Box const& bx, int ncomp,
                              Array4<Real> const& dUdt,
                              Array4<Real const> const& dUdt_in,
                              Array4<Real> const& scratch,
                              Array4<EBCellFlag const> const& flag,
                              Array4<Real const> const& vfrac)
{
    const Box dbox = Geom(lev).growPeriodicDomain(2);

    Array4<Real> tmp(scratch, 0);
    Array4<Real> delm(scratch, ncomp);
    Array4<Real> wgt(scratch, 2*ncomp);

    Box const& bxg1 = amrex::grow(bx,1);
    Box const& bxg2 = amrex::grow(bx,2);

    // xxxxx TODO: more weight options
    amrex::ParallelFor(bxg2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        wgt(i,j,k) = (dbox.contains(IntVect(i,j,k))) ? 1.0 : 0.0;
    });

    amrex::ParallelFor(bxg1, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isSingleValued()) {
            Real vtot = 0.0;
            Real divnc = 0.0;
            for (int kk = -1; kk <= 1; ++kk) {
            for (int jj = -1; jj <= 1; ++jj) {
            for (int ii = -1; ii <= 1; ++ii) {
                if ((ii != 0 or jj != 0 or kk != 0) and
                    flag(i,j,k).isConnected(ii,jj,kk) and
                    dbox.contains(IntVect(i+ii,j+jj,k+kk)))
                {
                    Real vf = vfrac(i+ii,j+jj,k+kk);
                    vtot += vf;
                    divnc += vf * dUdt_in(i+ii,j+jj,k+kk,n);
                }
            }}}
            divnc /= (vtot + 1.e-80);
            Real optmp = (1.0-vfrac(i,j,k))*(divnc-dUdt_in(i,j,k,n));
            tmp(i,j,k,n) = optmp;
            delm(i,j,k,n) = -vfrac(i,j,k)*optmp;
        } else {
            tmp(i,j,k,n) = 0.0;
        }
    });

    amrex::ParallelFor(bxg1 & dbox, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isSingleValued()) {
            Real wtot = 0.0;
            for (int kk = -1; kk <= 1; ++kk) {
            for (int jj = -1; jj <= 1; ++jj) {
            for (int ii = -1; ii <= 1; ++ii) {
                if ((ii != 0 or jj != 0 or kk != 0) and
                    flag(i,j,k).isConnected(ii,jj,kk))
                {
                    wtot += vfrac(i+ii,j+jj,k+kk) * wgt(i+ii,j+jj,k+kk);
                }
            }}}
            wtot = 1.0/(wtot+1.e-80);

            Real dtmp = delm(i,j,k,n) * wtot;
            for (int kk = -1; kk <= 1; ++kk) {
            for (int jj = -1; jj <= 1; ++jj) {
            for (int ii = -1; ii <= 1; ++ii) {
                if ((ii != 0 or jj != 0 or kk != 0) and
                    bx.contains(IntVect(i+ii,j+jj,k+kk)) and
                    flag(i,j,k).isConnected(ii,jj,kk))
                {
                    Gpu::Atomic::Add(&tmp(i+ii,j+jj,k+kk,n), dtmp*wgt(i+ii,j+jj,k+kk));
                }
            }}}
        }
    });

    amrex::ParallelFor(bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        dUdt(i,j,k,n) = dUdt_in(i,j,k,n) + tmp(i,j,k,n);
    });
}
#endif
