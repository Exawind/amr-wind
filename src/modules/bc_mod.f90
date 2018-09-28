!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: bc                                                     !
!  Purpose: Global variables for specifying boundary conditions.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module bc

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int

  use param, only: dim_bc

  ! Type of boundary:
  character(len=16) :: BC_Type(dim_bc)

  ! Flags for periodic boundary conditions
  logical :: cyclic_x = .false.
  logical :: cyclic_y = .false.
  logical :: cyclic_z = .false.

  ! Boundary condition coordinates
  real(rt) :: BC_X_w(dim_bc), BC_X_e(dim_bc)
  real(rt) :: BC_Y_s(dim_bc), BC_Y_n(dim_bc)
  real(rt) :: BC_Z_b(dim_bc), BC_Z_t(dim_bc)

  real(rt) :: BC_Normal(dim_bc,3)
  real(rt) :: BC_Center(dim_bc,3)

  ! Gas phase BC pressure
  real(rt) :: BC_P(dim_bc)

  ! Velocities at a specified boundary
  real(rt) :: BC_U(dim_bc)
  real(rt) :: BC_V(dim_bc)
  real(rt) :: BC_W(dim_bc)

  ! Volumetric flow rate through a mass inflow boundary
  real(rt) :: BC_VolFlow(dim_bc)

  ! Mass flow rate through a mass inflow boundary
  real(rt) :: BC_MassFlow(dim_bc)

  ! Specified pressure drop cyclic boundary
  real(rt) :: delp_x, delp_y, delp_z

  ! Partial slip wall boundary condition (gas only)
  real(rt) :: BC_hw(dim_bc)
  real(rt) :: BC_Uw(dim_bc)
  real(rt) :: BC_Vw(dim_bc)
  real(rt) :: BC_Ww(dim_bc)

  ! Heat transfer boundary condition
  real(rt) :: BC_T   (dim_bc)
  real(rt) :: BC_hw_T(dim_bc)
  real(rt) :: BC_Tw  (dim_bc)
  real(rt) :: BC_C_T (dim_bc)

  ! External shaking (shaking amplitude vector sets direction of shaking)
  real(rt), dimension(3) :: BC_shaker_A ! shaking amplitude
  real(rt)               :: BC_shaker_F ! shaking frequency

  ! Character variable to determine the flow plane of a flow cell
  character :: BC_Plane(dim_bc)

  ! Cell flag definitions
  integer, parameter :: undef_cell =   0 ! undefined
  integer, parameter :: pinf_      =  10 ! pressure inflow cell
  integer, parameter :: pout_      =  11 ! pressure outflow cell
  integer, parameter :: minf_      =  20 ! mass flux inflow cell
  integer, parameter :: nsw_       = 100 ! wall with no-slip b.c.
  integer, parameter :: fsw_       = 101 ! wall with free-slip
  integer, parameter :: psw_       = 102 ! wall with partial-slip b.c.
  integer, parameter :: cycl_      =  50 ! cyclic b.c.
  integer, parameter :: cycp_      =  51 ! cyclic b.c. with pressure drop


contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: set_cyclic                                               !
!                                                                      !
! Purpose: Function to set cyclic flags.                               !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  subroutine set_cyclic(cyc_x, cyc_y, cyc_z) &
    bind(C, name="incflo_set_cyclic")

    integer, intent(in) :: cyc_x, cyc_y, cyc_z

    cyclic_x = (cyc_x == 1)
    cyclic_y = (cyc_y == 1)
    cyclic_z = (cyc_z == 1)

  end subroutine set_cyclic


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: bc_defined                                               !
!                                                                      !
! Purpose: Return if a BC region has been defined based on coordinates !
! defined in the input deck.                                           !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  logical function bc_defined(icv)

    use param, only: is_defined

    integer, intent(in) :: icv

    bc_defined = is_defined(bc_x_w(icv)) .or. is_defined(bc_x_e(icv)) .or. &
                 is_defined(bc_y_s(icv)) .or. is_defined(bc_y_n(icv)) .or. &
                 is_defined(bc_z_b(icv)) .or. is_defined(bc_z_t(icv))

! An IC is defined for restart runs only if it is a 'PATCH'.
    if(bc_type(icv) == 'DUMMY') bc_defined = .false.

  end function bc_defined


end module bc
