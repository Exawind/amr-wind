module check_bc_outflow_module

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int
  use param,         only: one, undefined, zero, is_undefined, is_defined, equal
  use error_manager,  only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar, ival


  implicit none
  private

  public check_bc_p_outflow

contains

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  ! Subroutine: CHECK_BC_P_OUTFLOW                                       !
  ! Author: J.Musser                                    Date: 01-Mar-14  !
  !                                                                      !
  ! Purpose: Provided a detailed error message on bc                     !
  !                                                                      !
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  subroutine check_bc_p_outflow(BCV)

    use bc       , only: bc_p

    integer, intent(in) :: bcv

    call init_err_msg("CHECK_BC_P_OUTFLOW")

    if (is_undefined(bc_p(bcv))) then
       write(err_msg,1000) trim(ivar('BC_P',bcv))
       call flush_err_msg(abort=.true.)
    endif

1000 format('Error 1000: Required input not specified: ',A,/&
        'Please correct the input deck.')

    ! Clean up and return.
    call finl_err_msg

  end subroutine check_bc_p_outflow

end module check_bc_outflow_module
