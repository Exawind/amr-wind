!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: check_inputs                                            !
!  Purpose: read and verify input data, open files                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine check_inputs(dt) bind(C, name="check_inputs")

  use amrex_fort_module, only : rt => amrex_real

  use error_manager, only: init_error_manager
  use check_gas_prop_module, only: check_gas_properties

  real(rt)  , intent(in) :: dt

  ! Initialize the error manager. This call occurs after the incflo.dat
  ! is read so that message verbosity can be set and the .LOG file
  ! can be opened.
  call init_error_manager

  call check_gas_properties

end subroutine check_inputs
