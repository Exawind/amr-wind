module calc_cell_module

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   private

   public :: calc_cell_ic

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: CALC_CELL                                          !
!  Purpose: calculate the i, j or k cell index for the corresponding   !
!     x y or z reactor location. the index returned depends on which   !
!     half of the i, j or k cell that the x, y, or z position          !
!     intersects                                                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   pure integer function calc_cell(location, dx)

      use param, only : half, one

      implicit none

      ! the x, y or z location
      real(rt), intent(in) :: location

      ! the cell lengths along the corresponding axis (dx, dy or dz)
      real(rt), intent(in) :: dx

      calc_cell = floor(location/dx + half) - one

   end function calc_cell

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_cell_ic                                            !
!  Purpose: calculate the i, j or k cell index for IC regions.         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_cell_ic(dx, dy, dz, x_w, y_s, z_b, x_e, y_n, z_t, &
                           i_w, i_e, j_s, j_n, k_b, k_t)

      use param, only: equal

      implicit none

      real(rt), intent(in   ) :: dx, dy, dz
      real(rt), intent(in   ) :: x_w, y_s, z_b, x_e, y_n, z_t

      integer,      intent(  out) :: i_w, j_s, k_b, i_e, j_n, k_t

      i_w = calc_cell(x_w, dx) + 1
      i_e = calc_cell(x_e, dx)

      j_s = calc_cell(y_s, dy) + 1
      j_n = calc_cell(y_n, dy)

      k_b = calc_cell(z_b, dz) + 1
      k_t = calc_cell(z_t, dz)

   end subroutine calc_cell_ic

end module calc_cell_module
