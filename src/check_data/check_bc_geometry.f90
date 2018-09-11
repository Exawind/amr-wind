module check_bc_geometry_module

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   use error_manager, only: init_err_msg, finl_err_msg, flush_err_msg
   use error_manager, only: err_msg, ivar, ival

! Global Parameters:
!---------------------------------------------------------------------//
      use param, only: zero, undefined, undefined_i, undefined_c
      use param, only: is_undefined, is_defined

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_GEOMETRY                                        !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Determine if BCs are "DEFINED" and that they contain the    !
! minimum amount of geometry data.                                     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      subroutine check_bc_geometry

! Global Variables:
!---------------------------------------------------------------------//
! User specified BC
      use bc, only: bc_type, bc_defined
! User specified: BC geometry
      use bc, only: bc_x_e, bc_x_w
      use bc, only: bc_y_n, bc_y_s
      use bc, only: bc_z_t, bc_z_b

! Global Parameters:
!---------------------------------------------------------------------//
! The max number of BCs.
      use param, only: dim_bc

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! loop/variable indices
      integer :: BCV, I
! Error flag
      logical :: RECOGNIZED_BC_TYPE
! Total number of valid BC types
      integer, PARAMETER :: DIM_BCTYPE = 12
! Valid boundary condition types
      CHARACTER(LEN=16), DIMENSION(1:DIM_BCTYPE) ::VALID_BC_TYPE = (/&
           'MASS_INFLOW     ', 'MI              ',&
           'P_INFLOW        ', 'PI              ',&
           'P_OUTFLOW       ', 'PO              ',&
           'FREE_SLIP_WALL  ', 'FSW             ',&
           'NO_SLIP_WALL    ', 'NSW             ',&
           'PAR_SLIP_WALL   ', 'PSW             '/)
!......................................................................!

      call init_err_msg("CHECK_BC_GEOMETRY")

      do bcv = 1, dim_bc

         if(bc_type(bcv)/=undefined_c .and. bc_type(bcv)/='DUMMY')then
            recognized_bc_type = .false.
            do i = 1, dim_bctype
                valid_bc_type(i) = trim(valid_bc_type(i))
                if(valid_bc_type(i) == bc_type(bcv)) then
                   recognized_bc_type = .true.
                   exit
                endif
            enddo

            if(.not.recognized_bc_type) then
               write(err_msg, 1100) trim(ivar('BC_TYPE',bcv)), &
                  trim(bc_type(bcv)), valid_bc_type
               call flush_err_msg(abort=.true.)
            endif
         endif

 1100 format('Error 1100: Illegal entry: ',A,' = ',A,/'Valid entries:',&
         ' ',7(/5X,A,2x,A),/5X,A)

         if(bc_defined(bcv)) then

            if(is_undefined(bc_x_w(bcv))) then
               write(err_msg,1101) bcv, 'BC_X_w'
               call flush_err_msg(abort=.true.)
            endif

            if(is_undefined(bc_x_e(bcv))) then
               write(err_msg, 1101) bcv, 'BC_X_e'
               call flush_err_msg(abort=.true.)
            endif

            if(is_undefined(bc_y_s(bcv))) then
               write(err_msg, 1101) bcv, 'BC_Y_s'
               call flush_err_msg(abort=.true.)
            endif

            if(is_undefined(bc_y_n(bcv))) then
               write(err_msg, 1101) bcv, 'BC_Y_n'
               call flush_err_msg(abort=.true.)
            endif

            if(is_undefined(bc_z_b(bcv))) then
               write(err_msg, 1101) bcv, 'BC_Z_b'
               call flush_err_msg(abort=.true.)
            endif

            if(is_undefined(bc_z_t(bcv))) then
               write(err_msg, 1101) bcv, 'BC_Z_t'
               call flush_err_msg(abort=.true.)
            endif

 1101 FORMAT('Error 1101: Boundary condition ',I3,' is ill-defined.',/ &
         A,' are not specified.',/'Please correct the input deck.')

         endif
      enddo

      call finl_err_msg

      end subroutine check_bc_geometry

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: check_bc_geometry_wall                                  !
!                                                                      !
!  Purpose: Find and validate i, j, k locations for walls BC's         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine check_bc_geometry_wall(bcv, dx, dy, dz, &
         xlength, ylength, zlength, domlo, domhi)

! Global Variables:
!---------------------------------------------------------------------//
! Boundary condition locations and corresponding grid index
      use bc, only: bc_x_w, bc_x_e
      use bc, only: bc_y_s, bc_y_n
      use bc, only: bc_z_b, bc_z_t

! Function to compare two values
      use calc_cell_module, only: calc_cell_bc_wall

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//

      implicit none

      integer(c_int), intent(in) :: domlo(3),domhi(3)
      real(rt)  , intent(in) :: dx, dy, dz
      real(rt)  , intent(in) :: xlength, ylength, zlength

! Dummy Arguments:
!---------------------------------------------------------------------//
! Index of boundary condition.
      integer, INTENT(in) :: BCV

! Local Variables:
!---------------------------------------------------------------------//
! Calculated indices of the wall boundary
      integer :: I_w , I_e , J_s , J_n , K_b , K_t
! Integer error flag
      integer :: IER
!......................................................................!

      call init_err_msg("check_bc_geometry_wall")

      call calc_cell_bc_wall(domlo, domhi, &
         xlength, ylength, zlength, dx, dy, dz, &
         bc_x_w(bcv), bc_y_s(bcv), bc_z_b(bcv), &
         bc_x_e(bcv), bc_y_n(bcv), bc_z_t(bcv), &
         i_w, i_e, j_s, j_n, k_b, k_t)

! check for valid values
      ier = 0
      if (k_b<domlo(3)-1 .or. k_b>domhi(3)+1) ier = 1
      if (j_s<domlo(3)-1 .or. j_s>domhi(2)+1) ier = 1
      if (i_w<domlo(2)-1 .or. i_w>domhi(1)+1) ier = 1
      if (k_t<domlo(2)-1 .or. k_t>domhi(3)+1) ier = 1
      if (j_n<domlo(1)-1 .or. j_n>domhi(2)+1) ier = 1
      if (i_e<domlo(1)-1 .or. i_e>domhi(1)+1) ier = 1

      if (k_b > k_t) ier = 1
      if (j_s > j_n) ier = 1
      if (i_w > i_e) ier = 1

      if(ier /= 0)then
         write(err_msg,1100) bcv,                                      &
            'X', bc_x_w(bcv), bc_x_e(bcv), 'i', i_w, i_e, &
            'Y', bc_y_s(bcv), bc_y_n(bcv), 'j', j_s, j_n, &
            'Z', bc_z_b(bcv), bc_z_t(bcv), 'k', k_b, k_t
         call flush_err_msg(abort=.true.)
      endif

 1100 FORMAT('Error 1100: Invalid location specified for BC ',I3,'.',  &
         3(/3x,A1,': ',g12.5,',',g12.5,8x,A1,': ',I8,',',I8),/         &
         'Please correct the input deck.')

      call finl_err_msg

      end subroutine check_bc_geometry_wall

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: check_bc_geometry_flow                                  !
!                                                                      !
!  Purpose: Find and validate i, j, k locations for flow BC's. Also    !
!           set value of bc_plane for flow BC's.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine check_bc_geometry_flow(bcv, dx, dy, dz,&
         xlength, ylength, zlength, domlo, domhi)

! Global Variables:
!---------------------------------------------------------------------//
! Boundary condition locations and corresponding grid index
      use bc, only: BC_X_w, BC_X_e
      use bc, only: BC_Y_s, BC_Y_n
      use bc, only: BC_Z_b, BC_Z_t

      use calc_cell_module, only: calc_cell_bc_flow

      IMPLICIT NONE

      integer(c_int), intent(in) :: domlo(3),domhi(3)
      real(rt)  , intent(in) :: dx, dy, dz
      real(rt)  , intent(in) :: xlength, ylength, zlength

! Dummy Arguments:
!---------------------------------------------------------------------//
! Index of boundary condition.
      integer, INTENT(in) :: BCV

! Local Variables:
!---------------------------------------------------------------------//
! Calculated indices of the wall boundary
      integer :: I_w, I_e, J_s, J_n, K_b, K_t
! Indices for error checking
      integer :: ier

!......................................................................!

      call init_err_msg("CHECK_BC_GEOMETRY_FLOW")

      call calc_cell_bc_flow(&
         xlength, ylength, zlength, dx, dy, dz, &
         bc_x_w(bcv), bc_y_s(bcv), bc_z_b(bcv), &
         bc_x_e(bcv), bc_y_n(bcv), bc_z_t(bcv), &
         i_w, i_e, j_s, j_n, k_b, k_t)

! check for valid values
      ier = 0
      if(i_w<domlo(1)-1 .or. i_w>domhi(1)+1) ier = 1
      if(i_e<domlo(1)-1 .or. i_e>domhi(1)+1) ier = 1
      if(j_s<domlo(2)-1 .or. j_s>domhi(2)+1) ier = 1
      if(j_n<domlo(2)-1 .or. j_n>domhi(2)+1) ier = 1
      if(k_b<domlo(3)-1 .or. k_b>domhi(3)+1) ier = 1
      if(k_t<domlo(3)-1 .or. k_t>domhi(3)+1) ier = 1
      if(k_b > k_t) ier = 1
      if(j_s > j_n) ier = 1
      if(i_w > i_e) ier = 1

     if(ier /= 0)then
         write(err_msg,1100) bcv,                                      &
            'X', bc_x_w(bcv), bc_x_e(bcv), 'I', i_w, i_e, &
            'Y', bc_y_s(bcv), bc_y_n(bcv), 'J', j_s, j_n, &
            'Z', bc_z_b(bcv), bc_z_t(bcv), 'K', k_b, k_t
         call flush_err_msg(abort=.true.)
     endif

 1100 FORMAT('Error 1100: Invalid location specified for BC ',I3,'.',  &
         3(/3x,A1,': ',g12.5,',',g12.5,8x,A1,': ',I8,',',I8),/         &
         'Please correct the input deck.')

      call finl_err_msg

      end subroutine check_bc_geometry_flow

end module check_bc_geometry_module
