module eb_gradu_module

   use amrex_fort_module,  only: rt=>amrex_real, c_int
   use amrex_error_module, only: amrex_abort
   use param,              only: zero, half, one, two

   implicit none

   private
   public      :: compute_eb_gradu
   !
   ! We use no-slip boundary for velocities.
   !

contains

   subroutine compute_eb_gradu(gradu, dx, i, j, k, &
                                 vel, vlo, vhi,     &
                                 bcent, blo, bhi,   &
                                 apx, axlo, axhi,   &
                                 apy, aylo, ayhi,   &
                                 apz, azlo, azhi,   & 
                                 extrap_to_cc)

      ! Wall divergence operator
      real(rt),       intent(  out) :: gradu(9)

      ! Cell indeces
      integer(c_int), intent(in   ) :: i, j, k

      ! Grid spacing
      real(rt),       intent(in   ) :: dx(3)

      ! Array bounds
      integer(c_int), intent(in   ) ::  vlo(3),  vhi(3)
      integer(c_int), intent(in   ) :: axlo(3), axhi(3)
      integer(c_int), intent(in   ) :: aylo(3), ayhi(3)
      integer(c_int), intent(in   ) :: azlo(3), azhi(3)
      integer(c_int), intent(in   ) ::  blo(3),  bhi(3)

      ! Arrays
      real(rt),       intent(in   ) ::                               &
           &   vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3),     &
           & bcent(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),3),     &
           & apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)),   &
           & apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3)),   &
           & apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))

      ! If true, extrapolate gradu to cell centre instead of on EB wall 
      integer(c_int),  intent(in   ) :: extrap_to_cc

      ! Local variable
      real(rt)   :: dxinv(3)
      real(rt)   :: dapx, dapy, dapz
      real(rt)   :: apnorm, apnorminv, anrmx, anrmy, anrmz
      real(rt)   :: xit, yit, zit, s
      real(rt)   :: bct(3), d1, d2, ddinv
      real(rt)   :: cxm, cx0, cxp, cym, cy0, cyp, czm, cz0, czp
      real(rt)   :: u1, v1, w1, u2, v2, w2, dudn, dvdn, dwdn
      real(rt)   :: dist
      integer    :: ixit, iyit, izit, is

      ! Distance from wall to location where grad(u) is computed. 
      ! We leave it unless extrap_to_zero is true
      dist = zero

      dxinv = one / dx

      ! Difference in area fraction across cell, used to find 
      ! components of wall normal (AWAY from  wall), 
      ! and to compute the divergence at the end of subroutine
      dapx = apx(i+1,j,k)-apx(i,j,k)
      dapy = apy(i,j+1,k)-apy(i,j,k)
      dapz = apz(i,j,k+1)-apz(i,j,k)

      ! Magnitude of wall normal 
      apnorm = sqrt(dapx**2+dapy**2+dapz**2)

      if ( apnorm == zero ) then
         call amrex_abort("compute_diff_wallflux: we are in trouble.")
      end if

      ! Components of unit vector pointing TOWARD the wall 
      apnorminv = one/apnorm
      anrmx = -dapx * apnorminv
      anrmy = -dapy * apnorminv
      anrmz = -dapz * apnorminv

      ! The center of the wall, as an offset from the cell centre 
      bct = bcent(i,j,k,:)

      ! we want to move along the largest component of wall normal, 
      ! in order to interpolate from nearest adjacent cell centres
      if (abs(anrmx).ge.abs(anrmy) .and. abs(anrmx).ge.abs(anrmz)) then
         ! y-z plane: x = const
         ! the equation for the line:  x = bct(1) - d*anrmx
         !                             y = bct(2) - d*anrmy
         !                             z = bct(3) - d*anrmz
         if ( extrap_to_cc .eq. 1 ) then 
            dist = abs(bct(1))
         end if

         ! make sure we move in the right direction
         s = sign(one,-anrmx)
         is = nint(s)

         ! interpolation point 1: 
         ! the line intersects the y-z plane (x = s) at ...
         d1 = (bct(1) - s) * (one/anrmx)  ! this is also the distance from wall to intersection
         yit = bct(2) - d1*anrmy
         zit = bct(3) - d1*anrmz

         ! shift so that the center of the nine cells is located at (0.,0.)
         iyit = j + nint(yit)
         izit = k + nint(zit)
         yit = yit - nint(yit) 
         zit = zit - nint(zit)

         !
         ! coefficents for quadratic interpolation
         cym = half*yit*(yit-one)
         cy0 = one-yit*yit
         cyp = half*yit*(yit+one)
         czm = half*zit*(zit-one)
         cz0 = one-zit*zit
         czp = half*zit*(zit+one)

         !
         ! interpolation (note u,v,w = 0 at the wall)
         u1 = interp2d(cym,cy0,cyp,czm,cz0,czp, vel(i+is,iyit-1:iyit+1,izit-1:izit+1,1))
         v1 = interp2d(cym,cy0,cyp,czm,cz0,czp, vel(i+is,iyit-1:iyit+1,izit-1:izit+1,2))
         w1 = interp2d(cym,cy0,cyp,czm,cz0,czp, vel(i+is,iyit-1:iyit+1,izit-1:izit+1,3))

         !
         ! interpolation point 2: 
         ! the line intersects the y-z plane (x = 2*s) at ...
         d2 = (bct(1) - 2.d0*s) * (one/anrmx)
         yit = bct(2) - d2*anrmy
         zit = bct(3) - d2*anrmz

         ! shift so that the center of the nine cells are (0.,0.)
         iyit = j + nint(yit)
         izit = k + nint(zit)
         yit = yit - nint(yit)  
         zit = zit - nint(zit)

         !
         ! coefficents for quadratic interpolation
         cym = half*yit*(yit-one)
         cy0 = one-yit*yit
         cyp = half*yit*(yit+one)
         czm = half*zit*(zit-one)
         cz0 = one-zit*zit
         czp = half*zit*(zit+one)

         !
         ! interpolation
         u2 = interp2d(cym,cy0,cyp,czm,cz0,czp, vel(i+2*is,iyit-1:iyit+1,izit-1:izit+1,1))
         v2 = interp2d(cym,cy0,cyp,czm,cz0,czp, vel(i+2*is,iyit-1:iyit+1,izit-1:izit+1,2))
         w2 = interp2d(cym,cy0,cyp,czm,cz0,czp, vel(i+2*is,iyit-1:iyit+1,izit-1:izit+1,3))

      else if (abs(anrmy).ge.abs(anrmx) .and. abs(anrmy).ge.abs(anrmz)) then
         ! z-x plane

         if ( extrap_to_cc .eq. 1 ) then 
            dist = abs(bct(2))
         end if

         s = sign(one,-anrmy)
         is = nint(s)

         d1 = (bct(2) - s) * (one/anrmy)
         xit = bct(1) - d1*anrmx
         zit = bct(3) - d1*anrmz
         ixit = i + nint(xit)
         izit = k + nint(zit)
         xit = xit - nint(xit)
         zit = zit - nint(zit)

         cxm = half*xit*(xit-one)
         cx0 = one-xit*xit
         cxp = half*xit*(xit+one)
         czm = half*zit*(zit-one)
         cz0 = one-zit*zit
         czp = half*zit*(zit+one)

         u1 = interp2d(cxm,cx0,cxp,czm,cz0,czp, vel(ixit-1:ixit+1,j+is,izit-1:izit+1,1))
         v1 = interp2d(cxm,cx0,cxp,czm,cz0,czp, vel(ixit-1:ixit+1,j+is,izit-1:izit+1,2))
         w1 = interp2d(cxm,cx0,cxp,czm,cz0,czp, vel(ixit-1:ixit+1,j+is,izit-1:izit+1,3))

         d2 = (bct(2) - 2.d0*s) * (one/anrmy)
         xit = bct(1) - d2*anrmx
         zit = bct(3) - d2*anrmz
         ixit = i + nint(xit)
         izit = k + nint(zit)
         xit = xit - nint(xit)
         zit = zit - nint(zit)

         cxm = half*xit*(xit-one)
         cx0 = one-xit*xit
         cxp = half*xit*(xit+one)
         czm = half*zit*(zit-one)
         cz0 = one-zit*zit
         czp = half*zit*(zit+one)

         u2 = interp2d(cxm,cx0,cxp,czm,cz0,czp, vel(ixit-1:ixit+1,j+2*is,izit-1:izit+1,1))
         v2 = interp2d(cxm,cx0,cxp,czm,cz0,czp, vel(ixit-1:ixit+1,j+2*is,izit-1:izit+1,2))
         w2 = interp2d(cxm,cx0,cxp,czm,cz0,czp, vel(ixit-1:ixit+1,j+2*is,izit-1:izit+1,3))

      else
         ! x-y plane

         if ( extrap_to_cc .eq. 1 ) then 
            dist = abs(bct(3))
         end if

         s = sign(one,-anrmz)
         is = nint(s)

         d1 = (bct(3) - s) * (one/anrmz)
         xit = bct(1) - d1*anrmx
         yit = bct(2) - d1*anrmy
         ixit = i + nint(xit)
         iyit = j + nint(yit)
         xit = xit - nint(xit)
         yit = yit - nint(yit)

         cxm = half*xit*(xit-one)
         cx0 = one-xit*xit
         cxp = half*xit*(xit+one)
         cym = half*yit*(yit-one)
         cy0 = one-yit*yit
         cyp = half*yit*(yit+one)

         u1 = interp2d(cxm,cx0,cxp,cym,cy0,cyp, vel(ixit-1:ixit+1,iyit-1:iyit+1,k+is,1))
         v1 = interp2d(cxm,cx0,cxp,cym,cy0,cyp, vel(ixit-1:ixit+1,iyit-1:iyit+1,k+is,2))
         w1 = interp2d(cxm,cx0,cxp,cym,cy0,cyp, vel(ixit-1:ixit+1,iyit-1:iyit+1,k+is,3))

         d2 = (bct(3) - 2.d0*s) * (one/anrmz)
         xit = bct(1) - d2*anrmx
         yit = bct(2) - d2*anrmy
         ixit = i + nint(xit)
         iyit = j + nint(yit)
         xit = xit - nint(xit)
         yit = yit - nint(yit)

         cxm = half*xit*(xit-one)
         cx0 = one-xit*xit
         cxp = half*xit*(xit+one)
         cym = half*yit*(yit-one)
         cy0 = one-yit*yit
         cyp = half*yit*(yit+one)

         u2 = interp2d(cxm,cx0,cxp,cym,cy0,cyp, vel(ixit-1:ixit+1,iyit-1:iyit+1,k+2*is,1))
         v2 = interp2d(cxm,cx0,cxp,cym,cy0,cyp, vel(ixit-1:ixit+1,iyit-1:iyit+1,k+2*is,2))
         w2 = interp2d(cxm,cx0,cxp,cym,cy0,cyp, vel(ixit-1:ixit+1,iyit-1:iyit+1,k+2*is,3))

      end if

      !
      ! compute derivatives near the cell centre given that velocity is zero on the wall.
      !
      ddinv = one/(d1*d2*(d2-d1))
      dudn = -ddinv*(d2*d2*u1-d1*d1*u2 + two*(d1*u2-d2*u1)*dist)  ! note that the normal vector points toward the wall
      dvdn = -ddinv*(d2*d2*v1-d1*d1*v2 + two*(d1*v2-d2*v1)*dist)
      dwdn = -ddinv*(d2*d2*w1-d1*d1*w2 + two*(d1*w2-d2*w1)*dist)
      !
      ! transform them to d/dx, d/dy and d/dz given transverse derivatives are zero
      gradu(1) = dudn * anrmx
      gradu(2) = dudn * anrmy
      gradu(3) = dudn * anrmz
      !
      gradu(4) = dvdn * anrmx
      gradu(5) = dvdn * anrmy
      gradu(6) = dvdn * anrmz
      !
      gradu(7) = dwdn * anrmx
      gradu(8) = dwdn * anrmy
      gradu(9) = dwdn * anrmz

   end subroutine compute_eb_gradu

   real(rt) function interp2d(cym,cy0,cyp,czm,cz0,czp,v)
      real(rt), intent(in) :: cym,cy0,cyp,czm,cz0,czp,v(3,3)
      interp2d = czm*(cym*v(1,1) + cy0*v(2,1) + cyp*v(3,1)) &
                 +     cz0*(cym*v(1,2) + cy0*v(2,2) + cyp*v(3,2)) &
                 +     czp*(cym*v(1,3) + cy0*v(2,3) + cyp*v(3,3))
   end function interp2d

end module eb_gradu_module
