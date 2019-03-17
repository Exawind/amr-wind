module constant

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int, c_char

! Gravitational acceleration
   real(rt) :: gravity(3)

! Prescribed pressure gradient for cyclic BCs
   real(rt) :: delp(3)

! Initial density
   real(rt) :: ro_0

! Dynamic coefficient of viscosity
   real(rt) :: mu

! Flow index
   real(rt) :: n

! Yield stress
   real(rt) :: tau_0

! Papanastasiou regularisation parameter
   real(rt) :: papa_reg

! Zero-strain-limit viscosity
   real(rt) :: eta_0

! Fluid type
   character(:), allocatable :: fluid_model

! Initial conditions
   real(rt) :: ic_u
   real(rt) :: ic_v
   real(rt) :: ic_w
   real(rt) :: ic_p

! Parameters limiting maximum number of items for specifying BCs
   integer, parameter :: dim_bc = 6

! Parameters for testing if user input was specifed.
   real(rt),  parameter :: undefined = 9.87654321d31

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
      is_defined = not_equal(x, undefined)
   end function is_defined

   pure logical function not_equal(x, y)
      real(rt), intent(in) :: x, y
      not_equal = (abs(x-y) >= epsilon(x))
   end function not_equal

end module constant
