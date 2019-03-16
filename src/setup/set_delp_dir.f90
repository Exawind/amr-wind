subroutine set_delp_dir(delp_dir) &
   bind(C, name="set_delp_dir")

   use constant, only: delp, zero

   implicit none

   integer, intent(  out) :: delp_dir

   ! Pass out the direction in which we drop by delp (if any)
   ! so that we can set the correct periodicity flag in the C++
   if (abs(delp(1)) > epsilon(zero)) then
      delp_dir = 0
   else if (abs(delp(2)) > epsilon(zero)) then
      delp_dir = 1
   else if (abs(delp(3)) > epsilon(zero)) then
      delp_dir = 2
   else
      delp_dir = -1
   end if

end subroutine set_delp_dir
