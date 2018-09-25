!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module: fld_constants                                               C
!  Purpose: Common block containing field variable constants           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      module fld_const

         use amrex_fort_module, only : rt => amrex_real
         use iso_c_binding , only: c_int

         ! Specified constant gas density
         real(rt) :: ro_0

         ! Specified constant gas viscosity
         real(rt) :: mu_0

         ! Average molecular weight of gas
         real(rt) :: mw_avg

      end module fld_const
