!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: ic                                                          !
!  Purpose: Global initial conditions variables.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module ic

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   use param, only: dim_ic

   ! Boundary condition coordinates
   real(rt) :: IC_X_w(dim_ic), IC_X_e(dim_ic)
   real(rt) :: IC_Y_s(dim_ic), IC_Y_n(dim_ic)
   real(rt) :: IC_Z_b(dim_ic), IC_Z_t(dim_ic)

   ! Initial gas pressure
   real(rt) :: IC_P(dim_ic)

   ! Initial velocities in specified region
   real(rt) :: IC_U(dim_ic)
   real(rt) :: IC_V(dim_ic)
   real(rt) :: IC_W(dim_ic)

   ! Heat transfer boundary condition
   real(rt) :: IC_T(dim_ic)

   character(len=16) :: ic_pack_type(dim_ic)

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: ic_defined                                               !
!                                                                      !
! Purpose: Return if a IC region has been defined based on coordinates !
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
