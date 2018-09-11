module amrex_to_incflo_module
! _________________________________________________________________

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int, c_char
  implicit none
contains
!**************************************************************************!
!                                                                          !
!                                                                          !
!**************************************************************************!
  subroutine incflo_get_data( ro_g0_in, mu_g0_in, gravity_in ) &
    bind(C, name="incflo_get_data")

    use constant , only: gravity
    use fld_const, only: ro_g0, mu_g0
    use get_data_module, only: get_data
    use param, only: is_undefined

    implicit none

    real(rt),   intent(in ) :: ro_g0_in, mu_g0_in
    real(rt),   intent(in ) :: gravity_in(3)

    call get_data()

    mu_g0 = mu_g0_in
    ro_g0 = ro_g0_in
    gravity(:) = gravity_in(:)
                                          
  end subroutine incflo_get_data

end module amrex_to_incflo_module
