module calc_mu_module

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  subroutine: calc_mu
!                                                                      !
!  Purpose: Calculate the viscosity 
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine calc_mu(slo, shi, lo, hi, mu, lambda)

    use constant, only: mu_0

    use param, only: is_undefined

    implicit none

! Dummy arguments ....................................................//
    integer(c_int), intent(in   ) :: slo(3), shi(3), lo(3), hi(3)

    real(rt), intent(  out) ::  &
             mu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
         lambda(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Local variables .....................................................//
      integer :: i,j,k
      real(rt) :: mu_val, lambda_val

      ! Set the initial viscosity
      mu_val = mu_0
      lambda_val = -(2.0d0/3.0d0) * mu_val

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               mu(i,j,k) = mu_val
               lambda(i,j,k) = lambda_val
            enddo
         enddo
      enddo

    end subroutine calc_mu

  end module calc_mu_module
