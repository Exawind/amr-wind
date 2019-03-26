!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: bc                                                     !
!  Purpose: Global variables for specifying boundary conditions.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module bc

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   use constant, only: dim_bc, undefined, zero

   ! Type of boundary:
   character(len=16) :: bc_type(dim_bc)

   ! Flags for periodic boundary conditions
   logical :: cyclic_x = .false.
   logical :: cyclic_y = .false.
   logical :: cyclic_z = .false.

   logical :: bc_defined(1:dim_bc) = .false.

   ! Boundary condition location (EB planes)
   real(rt) :: bc_normal(1:dim_bc,1:3) = undefined
   real(rt) :: bc_center(1:dim_bc,1:3) = undefined

   ! Gas phase BC pressure
   real(rt) :: bc_p(1:dim_bc) = undefined

   ! Velocities at a specified boundary
   real(rt) :: bc_u(1:dim_bc) = zero
   real(rt) :: bc_v(1:dim_bc) = zero
   real(rt) :: bc_w(1:dim_bc) = zero

   ! Character variable to determine the flow plane of a flow cell
   character :: bc_plane(dim_bc)

   ! Cell flag definitions
   integer, parameter :: undef_cell =   0 ! undefined
   integer, parameter :: pinf_      =  10 ! pressure inflow cell
   integer, parameter :: pout_      =  11 ! pressure outflow cell
   integer, parameter :: minf_      =  20 ! mass flux inflow cell
   integer, parameter :: nsw_       = 100 ! wall with no-slip b.c.

contains

end module bc
