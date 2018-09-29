!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: incflo_get_walls                                          !
!                                                                      !
!  Author: J. Musser                                  Date: 05-FEB-17  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine incflo_get_walls(bcv, exists, normal, center) &
   bind(c,name='incflo_get_walls')

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   use bc, only: bc_defined, bc_type, bc_plane
   use bc, only: bc_normal, bc_center
   use bc, only: nsw_, fsw_, psw_

   use bc, only: bc_x_w, bc_y_s, bc_z_b
   use bc, only: bc_x_e, bc_y_n, bc_z_t
   use param, only: is_defined

   implicit none

   integer(c_int), intent(in   ) :: bcv
   integer(c_int), intent(  out) :: exists

   real(rt),   intent(  out) :: normal(3), center(3)

   real(rt) :: x, y, z

   real(rt), parameter :: offset = 1.0d-15

   exists = 0; 
   if (bc_defined(bcv)) then

      select case (trim(bc_type(bcv)))

      case('FREE_SLIP_WALL','FSW', &
           'NO_SLIP_WALL'  ,'NSW', &
           'PAR_SLIP_WALL' ,'PSW', &
           'MASS_INFLOW'   ,'MI')

         exists = 1; 
         ! Hack to override default plane orientation
         if(is_defined(bc_center(bcv,1))) then

            normal = bc_normal(bcv,:)
            center = bc_center(bcv,:)

         else
            select case (trim(bc_plane(bcv)))
            case('E'); 
               x = bc_x_w(bcv) + offset
               y = bc_y_s(bcv) + 0.5d0*(bc_y_s(bcv) + bc_y_n(bcv))
               z = bc_z_b(bcv) + 0.5d0*(bc_z_b(bcv) + bc_z_t(bcv))
               normal = (/ 1.0d0, 0.0d0, 0.0d0/)
            case('W'); 
               x = bc_x_e(bcv) - offset
               y = bc_y_s(bcv) + 0.5d0*(bc_y_s(bcv) + bc_y_n(bcv))
               z = bc_z_b(bcv) + 0.5d0*(bc_z_b(bcv) + bc_z_t(bcv))
               normal = (/-1.0d0, 0.0d0, 0.0d0/)
            case('N'); 
               x = bc_x_w(bcv) + 0.5d0*(bc_x_w(bcv) + bc_x_e(bcv))
               y = bc_y_s(bcv) + offset
               z = bc_z_b(bcv) + 0.5d0*(bc_z_b(bcv) + bc_z_t(bcv))
               normal = (/ 0.0d0, 1.0d0, 0.0d0/)
            case('S'); 
               x = bc_x_w(bcv) + 0.5d0*(bc_x_w(bcv) + bc_x_e(bcv))
               y = bc_y_n(bcv) - offset
               z = bc_z_b(bcv) + 0.5d0*(bc_z_b(bcv) + bc_z_t(bcv))
               normal = (/ 0.0d0,-1.0d0, 0.0d0/)
            case('T'); 
               x = bc_x_w(bcv) + 0.5d0*(bc_x_w(bcv) + bc_x_e(bcv))
               y = bc_y_s(bcv) + 0.5d0*(bc_y_s(bcv) + bc_y_n(bcv))
               z = bc_z_b(bcv) + offset
               normal = (/ 0.0d0, 0.0d0, 1.0d0/)
            case('B'); 
               x = bc_x_w(bcv) + 0.5d0*(bc_x_w(bcv) + bc_x_e(bcv))
               y = bc_y_s(bcv) + 0.5d0*(bc_y_s(bcv) + bc_y_n(bcv))
               z = bc_z_t(bcv) - offset
               normal = (/ 0.0d0, 0.0d0,-1.0d0/)
            end select
            center = (/ x, y, z/)

         endif

      end select
   endif

end subroutine incflo_get_walls

subroutine incflo_get_real_walls(bcv, exists, normal, center) &
   bind(c,name='incflo_get_real_walls')

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   use bc, only: bc_defined, bc_type, bc_plane
   use bc, only: bc_normal, bc_center
   use bc, only: nsw_, fsw_, psw_

   use bc, only: bc_x_w, bc_y_s, bc_z_b
   use bc, only: bc_x_e, bc_y_n, bc_z_t
   use param, only: is_defined

   implicit none

   integer(c_int), intent(in   ) :: bcv
   integer(c_int), intent(  out) :: exists

   real(rt),   intent(  out) :: normal(3), center(3)

   real(rt) :: x, y, z

   real(rt), parameter :: offset = 1.0d-15

   exists = 0; 
   if (bc_defined(bcv)) then

      select case (trim(bc_type(bcv)))

      case('FREE_SLIP_WALL','FSW', &
           'NO_SLIP_WALL'  ,'NSW', &
           'PAR_SLIP_WALL' ,'PSW')

         exists = 1; 
         ! Hack to override default plane orientation
         if(is_defined(bc_center(bcv,1))) then

            normal = bc_normal(bcv,:)
            center = bc_center(bcv,:)

         else
            select case (trim(bc_plane(bcv)))
            case('E'); 
               x = bc_x_w(bcv) + offset
               y = bc_y_s(bcv) + 0.5d0*(bc_y_s(bcv) + bc_y_n(bcv))
               z = bc_z_b(bcv) + 0.5d0*(bc_z_b(bcv) + bc_z_t(bcv))
               normal = (/ 1.0d0, 0.0d0, 0.0d0/)
            case('W'); 
               x = bc_x_e(bcv) - offset
               y = bc_y_s(bcv) + 0.5d0*(bc_y_s(bcv) + bc_y_n(bcv))
               z = bc_z_b(bcv) + 0.5d0*(bc_z_b(bcv) + bc_z_t(bcv))
               normal = (/-1.0d0, 0.0d0, 0.0d0/)
            case('N'); 
               x = bc_x_w(bcv) + 0.5d0*(bc_x_w(bcv) + bc_x_e(bcv))
               y = bc_y_s(bcv) + offset
               z = bc_z_b(bcv) + 0.5d0*(bc_z_b(bcv) + bc_z_t(bcv))
               normal = (/ 0.0d0, 1.0d0, 0.0d0/)
            case('S'); 
               x = bc_x_w(bcv) + 0.5d0*(bc_x_w(bcv) + bc_x_e(bcv))
               y = bc_y_n(bcv) - offset
               z = bc_z_b(bcv) + 0.5d0*(bc_z_b(bcv) + bc_z_t(bcv))
               normal = (/ 0.0d0,-1.0d0, 0.0d0/)
            case('T'); 
               x = bc_x_w(bcv) + 0.5d0*(bc_x_w(bcv) + bc_x_e(bcv))
               y = bc_y_s(bcv) + 0.5d0*(bc_y_s(bcv) + bc_y_n(bcv))
               z = bc_z_b(bcv) + offset
               normal = (/ 0.0d0, 0.0d0, 1.0d0/)
            case('B'); 
               x = bc_x_w(bcv) + 0.5d0*(bc_x_w(bcv) + bc_x_e(bcv))
               y = bc_y_s(bcv) + 0.5d0*(bc_y_s(bcv) + bc_y_n(bcv))
               z = bc_z_t(bcv) - offset
               normal = (/ 0.0d0, 0.0d0,-1.0d0/)
            end select
            center = (/ x, y, z/)

         endif

      end select
   endif

end subroutine incflo_get_real_walls
