!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: ic                                                          !
!  Purpose: Global initial conditions variables.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module ic

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   use param, only: dim_ic, undefined, zero

   ! Boundary condition coordinates
   real(rt) :: ic_x_w(1:dim_ic) = undefined
   real(rt) :: ic_x_e(1:dim_ic) = undefined
   real(rt) :: ic_y_s(1:dim_ic) = undefined
   real(rt) :: ic_y_n(1:dim_ic) = undefined
   real(rt) :: ic_z_b(1:dim_ic) = undefined
   real(rt) :: ic_z_t(1:dim_ic) = undefined

   ! Initial pressure
   real(rt) :: ic_p(1:dim_ic) = zero

   ! Initial velocities in specified region
   real(rt) :: ic_u(1:dim_ic) = zero
   real(rt) :: ic_v(1:dim_ic) = zero
   real(rt) :: ic_w(1:dim_ic) = zero

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: ic_defined                                               !
!                                                                      !
! Purpose: Return if a ic region has been defined based on coordinates !
! defined in the input deck.                                           !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   logical function ic_defined(icv)

      use param, only: is_defined

      integer, intent(in) :: icv

      ic_defined = &
         is_defined(ic_x_w(icv)) .or. is_defined(ic_x_e(icv)) .or. &
         is_defined(ic_y_s(icv)) .or. is_defined(ic_y_n(icv)) .or. &
         is_defined(ic_z_b(icv)) .or. is_defined(ic_z_t(icv))

   end function ic_defined

end module ic
