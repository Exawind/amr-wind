module check_boundary_conditions_module

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

! Parameter constants
   use param, only: zero, one, undefined, is_defined, is_undefined, equal

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
   use error_manager, only: init_err_msg, finl_err_msg, flush_err_msg
   use error_manager, only: err_msg, ivar, ival

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BOUNDARY_CONDITIONS                               !
!                                                                      !
!  Purpose: Check boundary condition specifications                    !
!     - convert physical locations to i, j, k's                        !
!     - convert mass and volumetric flows to velocities (FLOW_TO_VEL)  !
!     - check specification of physical quantities                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine check_boundary_conditions(dx, dy, dz, &
         xlength, ylength, zlength, domlo, domhi) &
         bind(C, name="check_boundary_conditions")

! Global Variables:
!---------------------------------------------------------------------//
! Total number of (actual) continuum solids.
      use constant, only: MMAX
! Flag: BC dimensions or Type is specified
      use bc, only: BC_DEFINED
! Use specified BC type
      use bc, only: BC_TYPE

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of BCs
      use param, only: DIM_BC
! Maximum number of disperse phases
      use param, only: DIM_M

      use check_bc_geometry_module, only: check_bc_geometry
      use check_bc_geometry_module, only: check_bc_geometry_flow
      use check_bc_geometry_module, only: check_bc_geometry_wall
      use check_bc_walls_module,    only: check_bc_walls
      use check_bc_inflow_module,   only: check_bc_p_inflow
      use check_bc_outflow_module,  only: check_bc_p_outflow

      implicit none

      integer(c_int), intent(in) :: domlo(3),domhi(3)
      real(rt)  , intent(in) :: dx,dy,dz
      real(rt)  , intent(in) :: xlength,ylength,zlength

! Local Variables:
!---------------------------------------------------------------------//
! Loop counters
      integer :: bcv
! Flag to skip checks on indexed solid phase.
      logical :: SKIP(1:DIM_M)
!......................................................................!

! Initialize the error manager.
      call init_err_msg("CHECK_BOUNDARY_CONDITIONS")

! Determine which BCs are DEFINED
      call check_bc_geometry

! Loop over each defined BC and check the user data.
      do bcv = 1, dim_bc

         if (bc_defined(bcv)) then

            select case (trim(bc_type(bcv)))

            case ('MASS_INFLOW', 'MI')
               call check_bc_geometry_flow(bcv,dx,dy,dz,&
                  xlength,ylength,zlength,domlo,domhi)

            case ('P_INFLOW', 'PI')
               call check_bc_geometry_flow(bcv,dx,dy,dz,&
                  xlength,ylength,zlength,domlo,domhi)
               call check_bc_p_inflow(mmax, skip, bcv)

            case ('P_OUTFLOW','PO')
               call check_bc_geometry_flow(bcv,dx,dy,dz,&
                  xlength,ylength,zlength,domlo,domhi)
               call check_bc_p_outflow(bcv)

            case ('FREE_SLIP_WALL','FSW',&
                  'PAR_SLIP_WALL', 'PSW',&
                  'NO_SLIP_WALL',  'NSW')
               call check_bc_geometry_wall(bcv,dx,dy,dz,&
                  xlength,ylength,zlength,domlo,domhi)
               call check_bc_walls(bcv)

            end select

         ! Check whether BC values are specified for undefined BC locations
         elseif(bc_type(bcv) /= 'DUMMY') then

            call check_bc_range(bcv)

         endif

      enddo

      ! Cleanup and exit.
      call finl_err_msg

   contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BC_RANGE                                          !
!                                                                      !
!  Purpose: Verify that data was not given for undefined BC regions.   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine check_bc_range(bcv)

      ! Gas phase BC varaibles
      use bc, only: BC_X_g, BC_P_g
      use bc, only: BC_U_g, BC_V_g, BC_W_g

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of disperse phases.
      use param, only: DIM_M
! Maximum number of species gas/solids
      use param, only: DIM_N_G, DIM_N_S

      use param, only: zero, one, equal

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------/
! Boundary condition index.
      integer, intent(in) :: bcv

! Local Variables:
!---------------------------------------------------------------------//
! Generic loop varaibles.
      integer :: N
!......................................................................!

! Initialize the error manager.
      call INIT_ERR_MSG("CHECK_BC_RANGE")

      if(.not.equal(bc_u_g(bcv),zero)) then
         write(err_msg,1100) trim(ivar('BC_U_g',bcv))
         call flush_err_msg(abort=.true.)
      endif

      if(.not.equal(bc_v_g(bcv),zero)) then
         write(err_msg,1100) trim(ivar('BC_V_g',bcv))
         call flush_err_msg(abort=.true.)
      endif

      if (.not.equal(bc_w_g(bcv),zero)) then
         write(err_msg,1100) trim(ivar('BC_W_g',bcv))
         call flush_err_msg(abort=.true.)
      endif

      if (is_defined(bc_p_g(bcv))) then
         write(err_msg,1100) trim(ivar('BC_P_g',bcv))
         call flush_err_msg(abort=.true.)
      endif

      DO N = 1, DIM_N_G
         IF(IS_DEFINED(BC_X_G(bcv,N))) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_X_g',bcv,N))
            call FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO

      call finl_err_msg

 1100 FORMAT('Error 1100:',A,' specified for an undefined BC location')

      end subroutine check_bc_range

   end subroutine check_boundary_conditions
end module check_boundary_conditions_module
