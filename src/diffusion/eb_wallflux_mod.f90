module eb_wallflux_mod

   use amrex_fort_module,  only: rt=>amrex_real, c_int
   use eb_gradu_module,    only: compute_eb_gradu
   use param,              only: zero, one, two
   use rheology_module,    only: viscosity

   implicit none

   private
   public      :: compute_diff_wallflux

contains

   !
   ! We use no-slip boundary for velocities.
   !
   subroutine compute_diff_wallflux (divw, dx, i, j, k, &
                                     vel, vlo, vhi,     &
                                     bcent, blo, bhi,   &
                                     apx, axlo, axhi,   &
                                     apy, aylo, ayhi,   &
                                     apz, azlo, azhi,   &
                                     do_explicit_diffusion)

      ! Wall divergence operator
      real(rt),       intent(  out) :: divw(3)

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
      real(rt),       intent(in   ) ::                                 &
           &   vel( vlo(1): vhi(1), vlo(2): vhi(2), vlo(3): vhi(3),3), &
           & bcent( blo(1): bhi(1), blo(2): bhi(2), blo(3): bhi(3),3), &
           &   apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)  ), &
           &   apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3)  ), &
           &   apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3)  )

      ! If true  then we include all the diffusive terms in this explicit result
      ! If false then we include all only the off-diagonal terms here -- we do this
      !     by computing the full tensor then subtracting the diagonal terms
      integer(c_int),  intent(in   ) :: do_explicit_diffusion

      ! Local variable
      real(rt)   :: idx, idy, idz
      real(rt)   :: dapx, dapy, dapz
      real(rt)   :: ux, uy, uz, vx, vy, vz, wx, wy, wz
      real(rt)   :: gradu(9)
      real(rt)   :: tauxx, tauyy, tauzz, tauxy, tauxz, tauyx, tauyz, tauzx, tauzy
      real(rt)   :: strain, visc

      divw  = zero
      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)

      ! Compute the velocity gradients on the EB wall 
      ! 
      call compute_eb_gradu(gradu, dx, i, j, k, & 
                            vel, vlo, vhi, &
                            bcent, blo, bhi, &
                            apx, axlo, axhi, & 
                            apy, aylo, ayhi, & 
                            apz, azlo, azhi, 0)

      ux = gradu(1) * idx
      uy = gradu(2) * idy
      uz = gradu(3) * idz
      !
      vx = gradu(4) * idx
      vy = gradu(5) * idy
      vz = gradu(6) * idz
      !
      wx = gradu(7) * idx
      wy = gradu(8) * idy
      wz = gradu(9) * idz

      strain = sqrt(two * ux**2 + two * vy**2 + two * wz**2 + (uy + vx)**2 + (vz + wy)**2 + (wx + uz)**2) 
      visc = viscosity(strain)

      ! compute components of stress tensor on the wall
      tauxx = visc * (ux + ux) 
      tauxy = visc * (uy + vx)
      tauxz = visc * (uz + wx)
      !
      tauyx = tauxy
      tauyy = visc * (vy + vy)
      tauyz = visc * (vz + wy)
      !
      tauzx = tauxz
      tauzy = tauyz
      tauzz = visc * (wz + wz)

      ! Difference in area fraction across cell, used to find the divergence
      dapx = apx(i+1,j,k)-apx(i,j,k)
      dapy = apy(i,j+1,k)-apy(i,j,k)
      dapz = apz(i,j,k+1)-apz(i,j,k)

      if (do_explicit_diffusion .eq. 0) then
         !
         ! Subtract diagonal terms of stress tensor, to be obtained through
         ! implicit solve instead.
         !
         tauxx = tauxx - visc * ux
         tauxy = tauxy - visc * uy
         tauxz = tauxz - visc * uz

         tauyx = tauyx - visc * vx
         tauyy = tauyy - visc * vy
         tauyz = tauyz - visc * vz

         tauzx = tauzx - visc * wx
         tauzy = tauzy - visc * wy
         tauzz = tauzz - visc * wz
      end if

      ! Return the divergence of the stress tensor 
      divw(1) = dapx*tauxx + dapy*tauyx + dapz*tauzx
      divw(2) = dapx*tauxy + dapy*tauyy + dapz*tauzy
      divw(3) = dapx*tauxz + dapy*tauyz + dapz*tauzz

   end subroutine compute_diff_wallflux

end module eb_wallflux_mod
