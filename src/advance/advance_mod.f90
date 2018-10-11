!
!
!  This module contains the subroutines used directly in advance()
!
!
module advance_mod

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

   !
   ! Add forcing (acceleration) terms to velocity
   !
   subroutine add_forcing(lo, hi, vel, ulo, uhi, dt) bind(C)

      use constant, only: gravity

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: ulo(3), uhi(3)

      ! Time step width
      real(ar),       intent(in   ) :: dt

      ! Arrays
      real(ar),       intent(inout) :: &
         vel(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)

      ! Local variables
      integer(c_int)                :: i, j, k , n

      do n = 1, 3
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)

                  vel(i,j,k,n) = vel(i,j,k,n) + dt * gravity(n)

               end do
            end do
         end do
      end do

   end subroutine add_forcing

end module advance_mod
