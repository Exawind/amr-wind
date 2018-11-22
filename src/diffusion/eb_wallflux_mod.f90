module eb_wallflux_mod

   use amrex_fort_module,  only: rt=>amrex_real, c_int
   use amrex_error_module, only: amrex_abort
   use param,              only: zero, half, one
   use eb_gradu_module,    only: compute_eb_gradu

   implicit none

   private
   public      :: compute_diff_wallflux

contains

   !
   ! We use no-slip boundary for velocities.
   !
   subroutine compute_diff_wallflux (divw, dx, i, j, k, &
                                     vel, vlo, vhi,     &
                                     eta, slo, shi, &
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
      integer(c_int), intent(in   ) ::  slo(3),  shi(3)
      integer(c_int), intent(in   ) :: axlo(3), axhi(3)
      integer(c_int), intent(in   ) :: aylo(3), ayhi(3)
      integer(c_int), intent(in   ) :: azlo(3), azhi(3)
      integer(c_int), intent(in   ) ::  blo(3),  bhi(3)

      ! Arrays
      real(rt),       intent(in   ) ::                               &
           &   vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3),     &
           &    eta(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),       &
           & bcent(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),3),     &
           & apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)),   &
           & apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3)),   &
           & apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))

      ! If true  then we include all the diffusive terms in this explicit result
      ! If false then we include all only the off-diagonal terms here -- we do this
      !     by computing the full tensor then subtracting the diagonal terms
      integer(c_int),  intent(in   ) :: do_explicit_diffusion

      ! Local variable
      real(rt)   :: dxinv(3)
      real(rt)   :: dapx, dapy, dapz
      real(rt)   :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
      real(rt)   :: gradu(9)
      real(rt)   :: tauxx, tauyy, tauzz, tauxy, tauxz, tauyx, tauyz, tauzx, tauzy

      divw  = zero
      dxinv = one / dx

      ! Compute the velocity gradients on the EB wall 
      ! 
      call compute_eb_gradu(gradu, dx, i, j, k, & 
                            vel, vlo, vhi, &
                            bcent, blo, bhi, &
                            apx, axlo, axhi, & 
                            apy, aylo, ayhi, & 
                            apz, azlo, azhi, 0)
      dudx = gradu(1)
      dudy = gradu(2)
      dudz = gradu(3)
      !
      dvdx = gradu(4)
      dvdy = gradu(5)
      dvdz = gradu(6)
      !
      dwdx = gradu(7)
      dwdy = gradu(8)
      dwdz = gradu(9)

      ! compute components of stress tensor on the wall
      tauxx = eta(i,j,k) * (dudx + dudx) 
      tauxy = eta(i,j,k) * (dudy + dvdx)
      tauxz = eta(i,j,k) * (dudz + dwdx)
      !
      tauyx = tauxy
      tauyy = eta(i,j,k) * (dvdy + dvdy)
      tauyz = eta(i,j,k) * (dvdz + dwdy)
      !
      tauzx = tauxz
      tauzy = tauyz
      tauzz = eta(i,j,k) * (dwdz + dwdz)

      ! Difference in area fraction across cell, used to find the divergence
      dapx = apx(i+1,j,k)-apx(i,j,k)
      dapy = apy(i,j+1,k)-apy(i,j,k)
      dapz = apz(i,j,k+1)-apz(i,j,k)

      if (do_explicit_diffusion .eq. 0) then
         !
         ! Subtract diagonal terms of stress tensor, to be obtained through
         ! implicit solve instead.
         !
         tauxx = tauxx - eta(i,j,k) * dudx
         tauxy = tauxy - eta(i,j,k) * dudy
         tauxz = tauxz - eta(i,j,k) * dudz

         tauyx = tauyx - eta(i,j,k) * dvdx
         tauyy = tauyy - eta(i,j,k) * dvdy
         tauyz = tauyz - eta(i,j,k) * dvdz

         tauzx = tauzx - eta(i,j,k) * dwdx
         tauzy = tauzy - eta(i,j,k) * dwdy
         tauzz = tauzz - eta(i,j,k) * dwdz
      end if

      ! Return the divergence of the stress tensor 
      divw(1) = dxinv(1) * (dapx*tauxx + dapy*tauyx + dapz*tauzx)
      divw(2) = dxinv(2) * (dapx*tauxy + dapy*tauyy + dapz*tauzy)
      divw(3) = dxinv(3) * (dapx*tauxz + dapy*tauyz + dapz*tauzz)

   end subroutine compute_diff_wallflux

end module eb_wallflux_mod
