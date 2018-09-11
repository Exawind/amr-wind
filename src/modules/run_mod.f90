module run

! Modules
!---------------------------------------------------------------------//

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

!---------------------------------------------------------------------//

      ! Stop-time of the run.
      real(rt) :: tstop

end module run
