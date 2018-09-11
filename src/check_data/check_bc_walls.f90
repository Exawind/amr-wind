module check_bc_walls_module

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int
  use error_manager,  only: finl_err_msg, flush_err_msg, init_err_msg, &
                          & ivar, ival, err_msg

  implicit none
  private

  public check_bc_walls


contains

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  ! Subroutine: CHECK_BC_WALLS                                           !
  ! Author: J.Musser                                    Date: 01-Mar-14  !
  !                                                                      !
  ! Purpose: Driver routine to call checks for WALL BCs.                 !
  !                                                                      !
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  subroutine check_bc_walls(bcv)

    integer, intent(in) :: bcv

    ! Initialize the error manager.
    call init_err_msg("CHECK_BC_WALLS")

    ! Input checks for gas phase.
    call check_bc_walls_gas(bcv)

    call finl_err_msg

  end subroutine check_bc_walls


  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  ! Subroutine: CHECK_BC_WALLS_GAS                                       !
  ! Author: J.Musser                                    Date: 01-Mar-14  !
  !                                                                      !
  ! Purpose: Check user-input for gas phase WALL BC parameters.          !
  !                                                                      !
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  subroutine check_bc_walls_gas(bcv)

    use bc,     only: bc_type, bc_uw_g, bc_vw_g, bc_ww_g
    use param, only: is_undefined


    integer, intent(in) :: bcv
    ! Initialize the error manger.
    call init_err_msg("CHECK_BC_WALLS_GAS")

    ! The wall velocities are not needed for no-slip or free-slip
    if(BC_TYPE(BCV) == 'PAR_SLIP_WALL' .or. BC_TYPE(BCV) == 'PSW') then
       if(is_undefined(bc_uw_g(bcv))) then
          write(err_msg,1000) trim(ivar('BC_Uw_g',bcv))
          call flush_err_msg(abort=.true.)
       elseif(is_undefined(bc_vw_g(bcv))) then
          write(err_msg,1000) trim(ivar('BC_Vw_g',bcv))
          call flush_err_msg(abort=.true.)
       elseif(is_undefined(bc_ww_g(bcv))) then
          write(err_msg,1000) trim(ivar('BC_Ww_g',bcv))
          call flush_err_msg(abort=.true.)
       endif
    endif

    ! Clear the error manager.
    call finl_err_msg

1000 format('Error 1000: Required input not specified: ',A,/&
       'Please correct the input deck.')

  end subroutine check_bc_walls_gas

end module check_bc_walls_module
