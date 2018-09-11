module calc_mu_g_module

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  subroutine: set_ro_g                                                !
!                                                                      !
!  Purpose: Initialize the gas density.                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine calc_mu_g(slo, shi, lo, hi, mu_g, lambda_g)

    use fld_const, only: mu_g0

    use param, only: is_undefined

    implicit none

! Dummy arguments ....................................................//
    integer(c_int), intent(in   ) :: slo(3), shi(3), lo(3), hi(3)

    real(rt), intent(  out) ::  &
             mu_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
         lambda_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Local variables .....................................................//
      integer :: i,j,k
      real(rt) :: mu_val, lambda_val

      ! Set the initial viscosity
      mu_val = mu_g0
      lambda_val = -(2.0d0/3.0d0) * mu_val

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               mu_g(i,j,k) = mu_val
               lambda_g(i,j,k) = lambda_val
            enddo
         enddo
      enddo

    end subroutine calc_mu_g

  end module calc_mu_g_module
