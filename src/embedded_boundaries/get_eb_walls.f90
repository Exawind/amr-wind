subroutine incflo_get_real_walls(bcv, exists, normal, center) &

   bind(c,name='incflo_get_real_walls')

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   use bc, only: bc_defined, bc_type
   use bc, only: bc_normal, bc_center
   use constant, only: is_defined

   implicit none

   integer(c_int), intent(in   ) :: bcv
   integer(c_int), intent(  out) :: exists

   real(rt),   intent(  out) :: normal(3), center(3)

   exists = 0; 

   if (bc_defined(bcv)) then

      select case (trim(bc_type(bcv)))

      case('NO_SLIP_WALL'  ,'NSW')

         exists = 1; 
         ! Hack to override default plane orientation
         if(is_defined(bc_center(bcv,1))) then

            normal = bc_normal(bcv,:)
            center = bc_center(bcv,:)

         endif

      end select
   endif

end subroutine incflo_get_real_walls
