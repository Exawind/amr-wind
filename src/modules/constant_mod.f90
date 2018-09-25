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

! Universal gas constant; (Pa.m3/kmol.K)
      real(rt), parameter :: gas_const = 8314.56D0

! Pi, the ubiquitous irrational number
      real(rt), PARAMETER :: pi = 4.d0*atan(1.d0)

! Solids phase constants
!-----------------------------------------------------------------------
! Number of solids phases
      integer :: mmax = 0

      END MODULE constant
