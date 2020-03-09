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
    constexpr Real small_vel = 1.e-10;

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

    // ****************************************************************************
    // Predict to x-faces
    // ****************************************************************************
    if ((extdir_ilo and domain_ilo >= ubx.smallEnd(0)-1) or
        (extdir_ihi and domain_ihi <= ubx.bigEnd(0)))
    {
        amrex::ParallelFor(Box(ubx),
        [u,vcc,flag,fcx,ccc,extdir_ilo,extdir_ihi,extdir_jlo,extdir_jhi,extdir_klo,extdir_khi,
         domain_ilo,domain_ihi,domain_jlo,domain_jhi,domain_klo,domain_khi]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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

               Real cc_umax = amrex::max(vcc(i,j,k,0), vcc(i-1,j,k,0));
               Real cc_umin = amrex::min(vcc(i,j,k,0), vcc(i-1,j,k,0));

               // Compute slopes of component "0" of vcc
               const auto& slopes_eb_hi = incflo_slopes_extdir_eb(i,j,k,0,vcc,ccc,flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real upls = vcc(i  ,j,k,0) - delta_x * slopes_eb_hi[0]
                                          + delta_y * slopes_eb_hi[1]
                                          + delta_z * slopes_eb_hi[2];

               upls = amrex::max(amrex::min(upls, cc_umax), cc_umin);

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

               umns = amrex::max(amrex::min(umns, cc_umax), cc_umin);

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
        [u,vcc,flag,fcx,ccc]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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

               Real cc_umax = amrex::max(vcc(i,j,k,0), vcc(i-1,j,k,0));
               Real cc_umin = amrex::min(vcc(i,j,k,0), vcc(i-1,j,k,0));

               // Compute slopes of component "0" of vcc
               const auto slopes_eb_hi = incflo_slopes_eb(i,j,k,0,vcc,ccc,flag);

               Real upls = vcc(i  ,j,k,0) - delta_x * slopes_eb_hi[0]
                                          + delta_y * slopes_eb_hi[1]
                                          + delta_z * slopes_eb_hi[2];

               upls = amrex::max(amrex::min(upls, cc_umax), cc_umin);

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

               umns = amrex::max(amrex::min(umns, cc_umax), cc_umin);

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
        [v,vcc,flag,fcy,ccc,extdir_ilo,extdir_ihi,extdir_jlo,extdir_jhi,extdir_klo,extdir_khi,
         domain_ilo,domain_ihi,domain_jlo,domain_jhi,domain_klo,domain_khi]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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

               Real cc_vmax = amrex::max(vcc(i,j,k,1), vcc(i,j-1,k,1));
               Real cc_vmin = amrex::min(vcc(i,j,k,1), vcc(i,j-1,k,1));

               // Compute slopes of component "1" of vcc
               const auto& slopes_eb_hi = incflo_slopes_extdir_eb(i,j,k,1,vcc,ccc,flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real vpls = vcc(i,j  ,k,1) + delta_x * slopes_eb_hi[0]
                                          - delta_y * slopes_eb_hi[1]
                                          + delta_z * slopes_eb_hi[2];

               vpls = amrex::max(amrex::min(vpls, cc_vmax), cc_vmin);

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

               vmns = amrex::max(amrex::min(vmns, cc_vmax), cc_vmin);

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
        [v,vcc,flag,fcy,ccc] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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

               Real cc_vmax = amrex::max(vcc(i,j,k,1), vcc(i,j-1,k,1));
               Real cc_vmin = amrex::min(vcc(i,j,k,1), vcc(i,j-1,k,1));

               // Compute slopes of component "1" of vcc
               const auto slopes_eb_hi = incflo_slopes_eb(i,j,k,1,vcc,ccc,flag);

               Real vpls = vcc(i,j  ,k,1) + delta_x * slopes_eb_hi[0]
                                          - delta_y * slopes_eb_hi[1]
                                          + delta_z * slopes_eb_hi[2];

               vpls = amrex::max(amrex::min(vpls, cc_vmax), cc_vmin);

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
                                          
               vmns = amrex::max(amrex::min(vmns, cc_vmax), cc_vmin);

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
        [w,vcc,flag,fcz,ccc,extdir_ilo,extdir_ihi,extdir_jlo,extdir_jhi,extdir_klo,extdir_khi,
         domain_ilo,domain_ihi,domain_jlo,domain_jhi,domain_klo,domain_khi]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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

               Real cc_wmax = amrex::max(vcc(i,j,k,2), vcc(i,j,k-1,2));
               Real cc_wmin = amrex::min(vcc(i,j,k,2), vcc(i,j,k-1,2));

               // Compute slopes of component "2" of vcc
               const auto& slopes_eb_hi = incflo_slopes_extdir_eb(i,j,k,2,vcc,ccc,flag,
                                          extdir_ilo, extdir_ihi, domain_ilo, domain_ihi,
                                          extdir_jlo, extdir_jhi, domain_jlo, domain_jhi,
                                          extdir_klo, extdir_khi, domain_klo, domain_khi);

               Real wpls = vcc(i,j,k  ,2) + delta_x * slopes_eb_hi[0]
                                          + delta_y * slopes_eb_hi[1]
                                          - delta_z * slopes_eb_hi[2];

               wpls = amrex::max(amrex::min(wpls, cc_wmax), cc_wmin);

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

               wmns = amrex::max(amrex::min(wmns, cc_wmax), cc_wmin);

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
        [w,vcc,flag,fcz,ccc] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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

               Real cc_wmax = amrex::max(vcc(i,j,k,2), vcc(i,j,k-1,2));
               Real cc_wmin = amrex::min(vcc(i,j,k,2), vcc(i,j,k-1,2));

               // Compute slopes of component "2" of vcc
               const auto slopes_eb_hi = incflo_slopes_eb(i,j,k,2,vcc,ccc,flag);

               Real wpls = vcc(i,j,k  ,2) + delta_x * slopes_eb_hi[0]
                                          + delta_y * slopes_eb_hi[1]
                                          - delta_z * slopes_eb_hi[2];

               wpls = amrex::max(amrex::min(wpls, cc_wmax), cc_wmin);

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

               wmns = amrex::max(amrex::min(wmns, cc_wmax), cc_wmin);

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
#endif
