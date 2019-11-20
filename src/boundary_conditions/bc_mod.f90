!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: bc                                                     !
!  Purpose: Global variables for specifying boundary conditions.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module bc

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int, c_char, c_null_char

   use constant, only: undefined, zero, one

   ! Type of boundary:
   character(len=16) :: bc_type(6)

   ! Flags for periodic boundary conditions
   logical :: cyclic_x = .false.
   logical :: cyclic_y = .false.
   logical :: cyclic_z = .false.

   logical :: bc_defined(1:6) = .false.

   ! Boundary condition location (EB planes)
   real(rt) :: bc_normal(1:6,1:3) = undefined
   real(rt) :: bc_center(1:6,1:3) = undefined

   ! Gas phase BC pressure
   real(rt) :: bc_p(1:6) = undefined

   ! Velocities at a specified boundary
   real(rt) :: bc_u(1:6) = zero
   real(rt) :: bc_v(1:6) = zero
   real(rt) :: bc_w(1:6) = zero

   ! Density at a specified boundary
   real(rt) :: bc_r(1:6) = one

   ! Tracer at a specified boundary
   real(rt) :: bc_t(1:6) = one

   ! Character variable to determine the flow plane of a flow cell
   character :: bc_plane(6)

   ! Cell flag definitions
   integer, parameter :: undef_cell =   0 ! undefined
   integer, parameter :: pinf_      =  10 ! pressure inflow cell
   integer, parameter :: pout_      =  11 ! pressure outflow cell
   integer, parameter :: minf_      =  20 ! mass flux inflow cell
   integer, parameter :: nsw_       = 100 ! wall with no-slip b.c.
   integer, parameter :: slip_      = 101 ! slip wall
   integer, parameter :: wall_model_ = 102 ! inhomogeneous neumann wall model

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutines: getters                                                 !
!                                                                      !
! Purpose: Getters for the boundary conditions values                  !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  integer(c_int) function get_bc_defined(pID) bind(C)
    integer(c_int), intent(in) :: pID
    if(bc_defined(pID)) then
      get_bc_defined = 1
    else
      get_bc_defined = 0
    endif
    return
  end function get_bc_defined

  subroutine get_bc_type(pID, c_string) bind(C)
    integer(c_int), intent(in) :: pID
    character(len=1, kind=c_char), intent(inout) :: c_string(16)
    integer :: N,I
    N = len_trim(BC_Type(pID))
    do I=1,N
      c_string(I) = BC_Type(pID)(I:I)
    enddo
    c_string(N+1) = c_null_char
  end subroutine get_bc_type

  real(rt) function get_bc_u(pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_bc_u = bc_u(pID)
    return
  end function get_bc_u

  real(rt) function get_bc_v(pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_bc_v = bc_v(pID)
    return
  end function get_bc_v

  real(rt) function get_bc_w(pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_bc_w = bc_w(pID)
    return
  end function get_bc_w

  real(rt) function get_bc_r(pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_bc_r = bc_r(pID)
    return
  end function get_bc_r

  real(rt) function get_bc_t(pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_bc_t = bc_t(pID)
    return
  end function get_bc_t

  real(rt) function get_bc_p(pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_bc_p = bc_p(pID)
    return
  end function get_bc_p
  
  integer(c_int) function get_minf() bind(C)
    get_minf = minf_
    return
  end function get_minf

  integer(c_int) function get_pinf() bind(C)
    get_pinf = pinf_
    return
  end function get_pinf

  integer(c_int) function get_pout() bind(C)
    get_pout = pout_
    return
  end function get_pout

  subroutine get_domain_bc (domain_bc_out) bind(C)
    integer(c_int), intent(out)  :: domain_bc_out(6)
    integer :: bcv
    ! Default is that we reflect particles off domain boundaries if not periodic
    domain_bc_out(1:6) = 1
    if (cyclic_x) domain_bc_out(1:2) = 0
    if (cyclic_y) domain_bc_out(3:4) = 0
    if (cyclic_z) domain_bc_out(5:6) = 0

    do bcv = 1,6
       select case (trim(bc_type(bcv)))
         case ('P_OUTFLOW','PO')
            domain_bc_out(bcv) = 0
       end select
    end do
  end subroutine get_domain_bc

end module bc
