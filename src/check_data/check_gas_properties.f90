module check_gas_prop_module

  use error_manager, only: init_err_msg, finl_err_msg, flush_err_msg, &
                           err_msg, ivar, ival

  use param, only: is_undefined, is_defined, zero

  implicit none
  private

  public check_gas_properties

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_GAS_PHASE                                        !
!  Purpose: Check the gas phase input section                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_gas_properties

    use fld_const, only: mu_0, ro_0, mw_avg

    call init_err_msg("CHECK_GAS_PHASE")

    ! We now require that ro_0 be defined and positive in incflo.dat
    if (is_defined(ro_0)) then
       if (ro_0 < zero) then ! incompressible
          write(err_msg, 1001) 'ro_0', ival(ro_0)
          call flush_err_msg(abort=.true.)
       endif
    else
       write(err_msg, 1001) 'ro_0', ival(ro_0)
       call flush_err_msg(abort=.true.)
    endif

    if (is_defined(mw_avg))then
       write(err_msg, 1100) 'ro_0 is specified'
       call flush_err_msg
    endif

    if (mu_0 < zero) then
       write(err_msg,1001) 'MU_0', ival(mu_0)
       call flush_err_msg(abort=.true.)
    endif

    call finl_err_msg

1100 format('Message 2000: MW_AVG is not needed when ',A,'.')
1001 format('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the input deck.')

  end subroutine check_gas_properties
end module check_gas_prop_module
