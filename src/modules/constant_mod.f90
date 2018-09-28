!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module: constant                                                    C
!  Author: M. Syamlal                                 Date: 5-FEB-92   C
!                                                                      C
!  Purpose: Common block containing physical constants and constants   C
!           used in the numerical technique                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      MODULE constant

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

! Modules
!---------------------------------------------------------------------//

! Gravitational acceleration
      real(rt) :: gravity(3)

! Initial density 
      real(rt) :: ro_0

! Dynamic coefficient of viscosity
      real(rt) :: mu_0

      END MODULE constant
