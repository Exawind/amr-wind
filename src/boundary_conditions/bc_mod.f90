!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: bc                                                     !
!  Purpose: Global variables for specifying boundary conditions.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module bc

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   use param, only: undefined
   use param, only: dim_bc

   ! Type of boundary:
   character(len=16) :: BC_Type(dim_bc)

   ! Flags for periodic boundary conditions
   logical :: cyclic_x = .false.
   logical :: cyclic_y = .false.
   logical :: cyclic_z = .false.

   logical :: bc_defined(1:dim_bc) = .false.

   ! Boundary condition location (EB planes)
   real(rt) :: BC_Normal(1:dim_bc,1:3) = undefined
   real(rt) :: BC_Center(1:dim_bc,1:3) = undefined

   ! Gas phase BC pressure
   real(rt) :: BC_P(dim_bc)

   ! Velocities at a specified boundary
   real(rt) :: BC_U(dim_bc)
   real(rt) :: BC_V(dim_bc)
   real(rt) :: BC_W(dim_bc)

   ! Specified pressure drop cyclic boundary
   real(rt) :: delp_x, delp_y, delp_z

   ! Character variable to determine the flow plane of a flow cell
   character :: BC_Plane(dim_bc)

   ! Cell flag definitions
   integer, parameter :: undef_cell =   0 ! undefined
   integer, parameter :: pinf_      =  10 ! pressure inflow cell
   integer, parameter :: pout_      =  11 ! pressure outflow cell
   integer, parameter :: minf_      =  20 ! mass flux inflow cell
   integer, parameter :: nsw_       = 100 ! wall with no-slip b.c.

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

end module bc
