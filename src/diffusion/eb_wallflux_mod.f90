module eb_wallflux_mod

   use amrex_fort_module,  only: rt=>amrex_real, c_int
   use amrex_error_module, only: amrex_abort
   use amrex_mlebabeclap_3d_module

   use constant,           only: zero, half, one, two
   use rheology_module,    only: viscosity

   implicit none

   private
   public      :: compute_diff_wallflux

contains

   !
   ! We use no-slip boundary for velocities.
   !
   subroutine compute_diff_wallflux(divw, dx, i, j, k, &
                                    vel,    vlo,  vhi, &
                                    bcent,  blo,  bhi, &
                                    flag,   flo,  fhi, &
                                    apx,   axlo, axhi, &
                                    apy,   aylo, ayhi, &
                                    apz,   azlo, azhi, &
                                    vfrac, vflo, vfhi, & 
                                    cyl_speed)

      ! Wall divergence operator
      real(rt),       intent(  out) :: divw(3)

      ! Cell indices
      integer(c_int), intent(in   ) :: i, j, k

      ! Grid spacing
      real(rt),       intent(in   ) :: dx(3)

      ! Array bounds
      integer(c_int), intent(in   ) ::  vlo(3),  vhi(3)
      integer(c_int), intent(in   ) :: axlo(3), axhi(3)
      integer(c_int), intent(in   ) :: aylo(3), ayhi(3)
      integer(c_int), intent(in   ) :: azlo(3), azhi(3)
      integer(c_int), intent(in   ) :: vflo(3), vfhi(3)
      integer(c_int), intent(in   ) ::  blo(3),  bhi(3)
      integer(c_int), intent(in   ) ::  flo(3),  fhi(3)

      ! Arrays
      real(rt),       intent(in   ) ::                                 &
           &   vel( vlo(1): vhi(1), vlo(2): vhi(2), vlo(3): vhi(3),3), &
           & bcent( blo(1): bhi(1), blo(2): bhi(2), blo(3): bhi(3),3), &
           &   apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)  ), &
           &   apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3)  ), &
           &   apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3)  ), & 
           & vfrac(vflo(1):vfhi(1),vflo(2):vfhi(2),vflo(3):vfhi(3)  )

      integer(c_int), intent(in   ) :: flag(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

      ! Rotating cylinder 
      real(rt),       intent(in   ) :: cyl_speed

      ! Local variable
      real(rt)   :: dxinv(3), idx, idy, idz
      real(rt)   :: dapx, dapy, dapz
      real(rt)   :: apnorm, apnorminv, anrmx, anrmy, anrmz
      real(rt)   :: dudn, dvdn, dwdn
      real(rt)   :: ux, uy, uz, vx, vy, vz, wx, wy, wz
      real(rt)   :: tauxx, tauxy, tauxz, tauyx, tauyy, tauyz, tauzx, tauzy, tauzz
      real(rt)   :: strain, visc
      real(rt)   :: ub, vb, wb, theta

      divw  = zero
      dxinv = one / dx

      idx = dxinv(1)
      idy = dxinv(2)
      idz = dxinv(3)

      dapx = apx(i+1,j  ,k  ) - apx(i,j,k)
      dapy = apy(i  ,j+1,k  ) - apy(i,j,k)
      dapz = apz(i  ,j  ,k+1) - apz(i,j,k)

      apnorm = sqrt(dapx**2 + dapy**2 + dapz**2)

      if ( apnorm == zero ) then
         call amrex_abort("compute_diff_wallflux: we are in trouble.")
      end if

      apnorminv = one/apnorm
      anrmx = -dapx * apnorminv  ! unit vector pointing toward the wall
      anrmy = -dapy * apnorminv
      anrmz = -dapz * apnorminv

      ! Value on wall 
      ub = 0.d0
      vb = 0.d0
      wb = 0.d0
      theta = 0.d0
      if (cyl_speed > 0.d0) then
         theta = atan2(-anrmy, -anrmx)
         ub =   cyl_speed * sin(theta)
         vb = - cyl_speed * cos(theta)
      endif


      call compute_dphidn_3d(dudn, dxinv, i, j, k, &
                             vel(:,:,:,1), vlo, vhi, &
                             flag, flo, fhi, &
                             bcent(i,j,k,:), ub,  &
                             anrmx, anrmy, anrmz, vfrac(i,j,k))
      call compute_dphidn_3d(dvdn, dxinv, i, j, k, &
                             vel(:,:,:,2), vlo, vhi, &
                             flag, flo, fhi, &
                             bcent(i,j,k,:), vb,  &
                             anrmx, anrmy, anrmz, vfrac(i,j,k))
      call compute_dphidn_3d(dwdn, dxinv, i, j, k, &
                             vel(:,:,:,3), vlo, vhi, &
                             flag, flo, fhi, &
                             bcent(i,j,k,:), wb,  &
                             anrmx, anrmy, anrmz, vfrac(i,j,k))

      !
      ! transform them to d/dx, d/dy and d/dz given transverse derivatives are zero
      ux = dudn * anrmx * idx
      uy = dudn * anrmy * idy
      uz = dudn * anrmz * idz
      !                       
      vx = dvdn * anrmx * idx
      vy = dvdn * anrmy * idy
      vz = dvdn * anrmz * idz
      !                       
      wx = dwdn * anrmx * idx
      wy = dwdn * anrmy * idy
      wz = dwdn * anrmz * idz

      strain = sqrt(two * ux**2 + two * vy**2 + two * wz**2 + &
                    (uy + vx)**2 + (vz + wy)**2 + (wx + uz)**2) 
      visc = viscosity(strain)

      ! compute components of EXPLICIT PART OF stress tensor on the wall
      tauxx = visc * ux 
      tauxy = visc * vx
      tauxz = visc * wx
      
      tauyx = visc * uy
      tauyy = visc * vy
      tauyz = visc * wy
      
      tauzx = visc * uz
      tauzy = visc * vz
      tauzz = visc * wz

      ! Difference in area fraction across cell, used to find the divergence
      dapx = apx(i+1,j,k)-apx(i,j,k)
      dapy = apy(i,j+1,k)-apy(i,j,k)
      dapz = apz(i,j,k+1)-apz(i,j,k)

      ! Return the divergence of the stress tensor 
      divw(1) = dapx * tauxx + dapy * tauxy + dapz * tauxz
      divw(2) = dapx * tauyx + dapy * tauyy + dapz * tauyz
      divw(3) = dapx * tauzx + dapy * tauzy + dapz * tauzz

   end subroutine compute_diff_wallflux

end module eb_wallflux_mod
