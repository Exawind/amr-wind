module amrex_to_incflo_module
! _________________________________________________________________

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int, c_char
  implicit none
contains
!**************************************************************************!
! Take constants from incflo and make them available in the Fortran module !
! "constant", so that they can be accessed from f90 functions.             !
!**************************************************************************!
  subroutine incflo_get_data(gravity_in, ro_0_in, mu_0_in) &
    bind(C, name="incflo_get_data")

    use constant , only: gravity, ro_0, mu_0
    use get_data_module, only: get_data
    use param, only: is_undefined

    implicit none

    real(rt),   intent(in ) :: gravity_in(3)
    real(rt),   intent(in ) :: ro_0_in, mu_0_in

    call get_data()

    gravity(:) = gravity_in(:)
    ro_0 = ro_0_in
    mu_0 = mu_0_in
                                          
  end subroutine incflo_get_data

end module amrex_to_incflo_module
