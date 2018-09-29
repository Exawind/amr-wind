!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: set_delp_dir                                            !
!                                                                      !
!  Purpose: Set delp_dir                                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine set_delp_dir(delp_dir) &
         bind(C, name="set_delp_dir")

         use bc   , only: delp_x, delp_y, delp_z
         use param, only: zero

         implicit none

         integer, intent(  out) :: delp_dir

         ! Pass out the direction in which we drop by delp (if any)
         ! so that we can set the correct periodicity flag in the C++
         if (abs(delp_x) > epsilon(zero)) then
            delp_dir = 0
         else if (abs(delp_y) > epsilon(zero)) then
            delp_dir = 1
         else if (abs(delp_z) > epsilon(zero)) then
            delp_dir = 2
         else
            delp_dir = -1
         end if

      end subroutine set_delp_dir
