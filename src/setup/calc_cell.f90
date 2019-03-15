module calc_cell_module

   use amrex_fort_module, only : rt => amrex_real
   use param, only: half, one

   private

   public :: calc_cell_ic

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_cell_ic                                            !
!  Purpose: calculate the i, j or k cell index for IC regions.         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_cell_ic(dx, dy, dz, x_w, y_s, z_b, x_e, y_n, z_t, &
                           i_w, i_e, j_s, j_n, k_b, k_t)

      implicit none

      real(rt), intent(in   ) :: dx, dy, dz
      real(rt), intent(in   ) :: x_w, y_s, z_b, x_e, y_n, z_t

      integer,  intent(  out) :: i_w, j_s, k_b, i_e, j_n, k_t

      i_w = floor(x_w / dx + half)
      i_e = floor(x_e / dx + half) - one
      j_s = floor(y_s / dy + half)
      j_n = floor(y_n / dy + half) - one
      k_b = floor(z_b / dz + half)
      k_t = floor(z_t / dz + half) - one

   end subroutine calc_cell_ic

end module calc_cell_module
