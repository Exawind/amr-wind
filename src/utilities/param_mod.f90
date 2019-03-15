module param

   use amrex_fort_module, only : rt => amrex_real

! Parameters limiting maximum number of items for specifying ICs and BCs
   integer, parameter :: dim_ic = 6
   integer, parameter :: dim_bc = 6

! Parameters for testing if user input was specifed.
   real(rt),  parameter :: undefined = 9.87654321d31
   character, parameter :: undefined_c = ' '

! Cutoffs for large and small numbers
   real(rt), parameter :: my_huge  = 1.0d20

! Common parameter constants
   real(rt), parameter :: zero = 0.0d0
   real(rt), parameter :: half = 0.5d0
   real(rt), parameter :: one  = 1.0d0
   real(rt), parameter :: two  = 2.0d0

contains

   pure logical function is_defined(x)
      real(rt), intent(in) :: x
      is_defined = .not.equal(x, undefined)
   end function is_defined

   pure logical function is_undefined(x)
      real(rt), intent(in) :: x
      is_undefined = equal(x, undefined)
   end function is_undefined

   pure logical function equal(x, y)
      real(rt), intent(in) :: x, y
      equal = (abs(x-y) < epsilon(x))
   end function equal

end module param
