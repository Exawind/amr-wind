!
!
!  This module contains the subroutines to perform some of the steps of the
!  MAC projection method.
!
!
module macproj_mod

   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one

   implicit none
   private

   ! Define here the unit vectors
   ! This is used to shift index  based on how the variable is staggered
   ! Check e_x, e_y and e_z in incflo.H
   integer(c_int), parameter :: e_i(3,3) = reshape ( [1,0,0,0,1,0,0,0,1], [3,3] )

contains

   subroutine compute_bcoeff_mac ( lo, hi, bcoeff, blo, bhi, &
                                  ro, slo, shi, dir )  bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: blo(3),bhi(3)

      ! Direction
      integer(c_int), intent(in   ) :: dir

      ! Arrays
      real(ar),       intent(in   ) :: &
         ro(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),       intent(  out) :: &
         bcoeff(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3))

      integer      :: i, j, k, i0, j0, k0
      real(ar)     :: ro_f

      i0 = e_i(dir,1)
      j0 = e_i(dir,2)
      k0 = e_i(dir,3)

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               ro_f = half * ( ro(i,j,k) + ro(i-i0,j-j0,k-k0) )

               bcoeff(i,j,k) = one / ro_f

            end do
         end do
      end do

   end subroutine compute_bcoeff_mac

   !
   ! Computes  u_i = u_i + C * (1/ro) * (dphi/dx_i)
   !
   ! u_i      = i-th component of a staggered vector field u.
   !
   ! ro       = cell centered density field
   !
   ! dphidxi  =  i-th component of staggered gradphi
   !
   ! C        = real constant
   !
   ! dir      = 1,2,3 indicates x, y, z direction respectively.
   !
   subroutine project_mac_velocity ( lo, hi, u_i, ulo, uhi, &
        & dphidxi, glo, ghi, ro, slo, shi, c, dir ) bind (C)

      ! Loop bounds
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array bounds
      integer(c_int),  intent(in   ) :: ulo(3), uhi(3)
      integer(c_int),  intent(in   ) :: glo(3), ghi(3)
      integer(c_int),  intent(in   ) :: slo(3), shi(3)

      ! Grid and time spacing
      real(ar),        intent(in   ) :: c

      ! Direction
      integer(c_int),  intent(in   ) :: dir

      ! Arrays
      real(ar),        intent(in   ) ::                          &
           &      ro(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           & dphidxi(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3))

      real(ar),        intent(inout) ::                           &
           & u_i(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      ! Local variables
      integer(c_int)                 :: i, j, k, i0, j0, k0
      real(ar)                       :: oro

      i0 = e_i(dir,1)
      j0 = e_i(dir,2)
      k0 = e_i(dir,3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               oro  = half * ( one/ro(i,j,k) + one/ro(i-i0,j-j0,k-k0) )
               u_i(i,j,k) = u_i(i,j,k) + c * oro * dphidxi(i,j,k)
            end do
         end do
      end do

   end subroutine project_mac_velocity

end module macproj_mod
