!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: ic                                                          !
!  Purpose: Global initial conditions variables.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module ic

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   use param, only: undefined, zero

   ! Boundary condition coordinates
   real(rt) :: ic_x_w = undefined
   real(rt) :: ic_x_e = undefined
   real(rt) :: ic_y_s = undefined
   real(rt) :: ic_y_n = undefined
   real(rt) :: ic_z_b = undefined
   real(rt) :: ic_z_t = undefined

   ! Initial pressure
   real(rt) :: ic_p = zero

   ! Initial velocities in specified region
   real(rt) :: ic_u = zero
   real(rt) :: ic_v = zero
   real(rt) :: ic_w = zero

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: ic_defined                                               !
!                                                                      !
! Purpose: Return if a ic region has been defined based on coordinates !
! defined in the input deck.                                           !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   logical function ic_defined()

      use param, only: is_defined

      ic_defined = &
         is_defined(ic_x_w) .or. is_defined(ic_x_e) .or. &
         is_defined(ic_y_s) .or. is_defined(ic_y_n) .or. &
         is_defined(ic_z_b) .or. is_defined(ic_z_t)

   end function ic_defined

end module ic
