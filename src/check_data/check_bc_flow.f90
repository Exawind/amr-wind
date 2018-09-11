module check_bc_flow_module

   use param, only: zero
   use param,  only: dim_bc, dim_m

   use bc, only: bc_plane
   use bc, only: bc_u_g, bc_v_g, bc_w_g
   use bc, only: bc_u_s, bc_v_s, bc_w_s

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
   use error_manager, only: init_err_msg, finl_err_msg, flush_err_msg
   use error_manager, only: err_msg, ivar, ival

   implicit none

   private
   public check_bc_flow
contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: check_bc_flow                                           !
!                                                                      !
!  Purpose: Check boundary condition specifications                    !
!     - convert physical locations to i, j, k's                        !
!     - convert mass and volumetric flows to velocities (FLOW_TO_VEL)  !
!     - check specification of physical quantities                     !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine check_bc_flow() &
      bind(C,name ="check_bc_flow")

      use bc, only: bc_defined, bc_type

      implicit none

      integer :: bcv

      ! Loop over each defined BC and check the user data.
      do bcv = 1, dim_bc

         if (bc_defined(bcv)) then

            select case (trim(bc_type(bcv)))
            case ('MASS_INFLOW','MI')
               call check_bc_vel_inflow(bcv)

            case ('MASS_OUTFLOW','MO')
               call check_bc_vel_outflow(bcv)
            end select
         endif
      enddo

   end subroutine check_bc_flow

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: check_bc_vel_inflow                                      !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_bc_vel_inflow(bcv)

    use param, only: dim_m

    integer, intent(in) :: bcv

    ! Define format for error messages

    call init_err_msg("CHECK_BC_VEL_INFLOW")

    ! Check that gas phase velocities are consistent.
    select case (bc_plane(bcv))

    case ('W')
       if(bc_u_g(bcv) > zero) then
          write(err_msg,1300) trim(ivar('BC_U_g',bcv)), '<'
          call flush_err_msg(abort=.true.)
       endif
    case('E')
       if(bc_u_g(bcv) < zero) then
          write(err_msg,1300) trim(ivar('BC_U_g',bcv)), '>'
          call flush_err_msg
       endif
    case('S')
       if(bc_v_g(bcv) > zero) then
          write(err_msg,1300) trim(ivar('BC_V_g',bcv)), '<'
          call flush_err_msg
       endif
    case('N')
       if(bc_v_g(bcv) < zero) then
          write(err_msg,1300) trim(ivar('BC_V_g',bcv)), '>'
          call flush_err_msg
       endif
    case('B')
       if(bc_w_g(bcv) > zero) then
          write(err_msg,1300) trim(ivar('BC_W_g',bcv)), '<'
          call flush_err_msg
       endif
    case('T')
       if(bc_w_g(bcv) < zero) then
          write(err_msg,1300) trim(ivar('BC_W_g',bcv)), '>'
          call flush_err_msg
       endif

    end select

    call finl_err_msg

 1300 format('Error 1300: Invalid flow direction. '/,A,' should be ', &
          A,' zero. ',/'Please correct the input deck.')

  end subroutine check_bc_vel_inflow


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: check_bc_vel_outflow                                     !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_bc_vel_outflow(bcv)

    use bc, only: dim_m

    ! loop/variable indices
    integer, intent(in) :: bcv

    call init_err_msg("CHECK_BC_VEL_OUTFLOW")

    ! Check that gas phase velocities are consistent.
    select case (bc_plane(BCV))

    case ('W')
       if(bc_u_g(bcv) < zero) then
          write(err_msg,1300) trim(ivar('BC_U_g',bcv)), '>'
          call flush_err_msg
       endif
    case('E')
       if(bc_u_g(bcv) > zero) then
          write(err_msg,1300) trim(ivar('BC_U_g',bcv)), '<'
          call flush_err_msg
       endif
    case('S')
       if(bc_v_g(bcv) < zero) then
          write(err_msg,1300) trim(ivar('BC_V_g',bcv)), '>'
          call flush_err_msg
       endif
    case('N')
       if(bc_v_g(bcv) > zero) then
          write(err_msg,1300) trim(ivar('BC_V_g',bcv)), '<'
          call flush_err_msg
       endif
    case('B')
       if(bc_w_g(bcv) < zero) then
          write(err_msg,1300) trim(ivar('BC_W_g',bcv)), '>'
          call flush_err_msg
       endif
    case('T')
       if(bc_w_g(bcv) > zero) then
          write(err_msg,1300) trim(ivar('BC_W_g',bcv)), '<'
          call flush_err_msg
       endif

    end select

    call finl_err_msg

 1300 format('Error 1300: Invalid flow direction. ',/A,' should be ', &
          A,' zero. ',/'Please correct the input deck.')

  end subroutine check_bc_vel_outflow
end module check_bc_flow_module
