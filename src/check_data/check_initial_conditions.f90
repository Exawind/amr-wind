module check_initial_conditions_module

  use amrex_fort_module, only: rt => amrex_real
  use iso_c_binding ,    only: c_int

  use param,  only: zero, one, undefined, undefined_i
  use param,  only: is_defined, is_undefined

  use error_manager, only: init_err_msg, finl_err_msg, flush_err_msg
  use error_manager, only: err_msg, ivar, ival

  implicit none
  private

  public check_initial_conditions

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_INITIAL_CONDITIONS                                !
!                                                                      !
!  Purpose: check the initial conditions input section                 !
!     - check geometry of any specified IC region                      !
!     - check specification of physical quantities                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_initial_conditions(dx,dy,dz,domlo,domhi) &
      bind(C, name="check_initial_conditions")

    use ic,                    only: ic_defined
    use param,                 only: dim_ic

    integer(c_int), intent(in) :: domlo(3),domhi(3)
    real(rt)  , intent(in) :: dx, dy, dz
    integer(c_int)             :: icv

    ! Determine which ICs are DEFINED
    call check_ic_geometry(dx,dy,dz,domlo,domhi)

    ! Loop over all IC arrays.
    do icv=1, dim_ic

       ! Verify user input for defined defined IC.
       if(ic_defined(icv)) then
          ! Gas phase checks.
          call check_ic_gas_phase(ICV)

       ! Verify that no data was defined for unspecified IC.
       else
          call check_ic_overflow(ICV)
       endif
    enddo

 end subroutine check_initial_conditions

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_IC_GEOMETRY                                        !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    subroutine check_ic_geometry(dx,dy,dz,domlo,domhi)

      use param, only: dim_ic
      use ic,    only: ic_defined
      use ic,    only: IC_X_e, IC_Y_n, IC_Z_t
      use ic,    only: IC_X_w, IC_Y_s, IC_Z_b

      use calc_cell_module, only: calc_cell_ic

      integer(c_int), intent(in) :: domlo(3),domhi(3)
      real(rt)  , intent(in) :: dx,dy,dz

      integer :: icv, i_w, i_e, j_s, j_n, k_b, k_t

      ! Initialize the error manager.
      call init_err_msg("CHECK_IC_GEOMETRY")

      ! Check geometry of any specified IC region
      do icv = 1, dim_ic

         if(ic_defined(icv)) then

            if (is_undefined(ic_x_w(icv))) then
               write(err_msg, 1100) icv, 'IC_X_w'
               call flush_err_msg(abort=.true.)
            endif

           if (is_undefined(ic_x_e(icv))) then
              write(err_msg, 1100) icv, 'IC_X_e'
              call flush_err_msg(abort=.true.)
           endif

           if (is_undefined(ic_y_s(icv))) then
              write(err_msg, 1100) icv, 'IC_Y_s'
              call flush_err_msg(abort=.true.)
           endif

           if (is_undefined(ic_y_n(icv))) then
              write(err_msg, 1100) icv, 'IC_Y_n'
              call flush_err_msg(abort=.true.)
           endif

           if (is_undefined(ic_z_b(icv))) then
              write(err_msg, 1100) icv, 'IC_Z_b'
              call flush_err_msg(abort=.true.)
           endif

           if (is_undefined(ic_z_t(icv))) then
              write(err_msg, 1100) icv, 'IC_Z_t'
              call flush_err_msg(abort=.true.)
           endif

1100 FORMAT('Error 1100: Initial condition region ',I3,' is ill-',    &
        'defined.',/' > ',A,' is not specified.',/'Please correct ', &
        'the input deck.')

           call calc_cell_ic(dx, dy, dz, &
                  ic_x_w(icv), ic_y_s(icv), ic_z_b(icv), &
                  ic_x_e(icv), ic_y_n(icv), ic_z_t(icv), &
                  i_w, i_e, j_s, j_n, k_b, k_t)

           ! Report problems with calculated bounds.
           if(i_w > i_e) then
              write(err_msg, 1101) icv, 'IC_I_W > IC_I_E'
              call flush_err_msg(abort=.true.)
           elseif(i_w < domlo(1)) then
              write(err_msg, 1101) icv, 'IC_I_W < domlo(1)'
              call flush_err_msg(abort=.true.)
           elseif(i_w > domhi(1)) then
              write(err_msg, 1101) icv, 'IC_I_W > domhi(1)'
              call flush_err_msg(abort=.true.)
           elseif(i_e < domlo(1)) then
              write(err_msg, 1101) icv, 'IC_I_E < domlo(1)'
              call flush_err_msg(abort=.true.)
           elseif(i_e > domhi(1)) then
              write(err_msg, 1101) icv, 'IC_I_E > domhi(1)'
              call flush_err_msg(abort=.true.)
           endif

           if(j_s > j_n) then
              write(err_msg, 1101) icv, 'IC_J_S > IC_J_N'
              call flush_err_msg(abort=.true.)
           elseif(j_s<domlo(2)) then
              write(err_msg, 1101) icv, 'IC_J_S < domlo(2)'
              call flush_err_msg(abort=.true.)
           elseif(j_s>domhi(2)) then
              write(err_msg, 1101) icv, 'IC_J_S >  domhi(2)'
              call flush_err_msg(abort=.true.)
           elseif(j_n<domlo(2)) then
              write(err_msg, 1101) icv, 'IC_J_N < domlo(2)'
              call flush_err_msg(abort=.true.)
           elseif(j_n>domhi(2)) then
              write(err_msg, 1101) icv, 'IC_J_N > domhi(2)'
              call flush_err_msg(abort=.true.)
           endif

           if(k_b > k_t) then
              write(err_msg, 1101) icv, 'IC_K_B > IC_K_T'
              call flush_err_msg(abort=.true.)
           elseif(k_b < domlo(3)) then
              write(err_msg, 1101) icv, 'IC_K_B < domlo(3)'
              call flush_err_msg(abort=.true.)
           elseif(k_b > domhi(3)) then
              write(err_msg, 1101) icv, 'IC_K_B > domhi(3)'
              call flush_err_msg(abort=.true.)
           elseif(k_t < domlo(3)) then
              write(err_msg, 1101) icv, 'IC_K_T < domlo(3)'
              call flush_err_msg(abort=.true.)
           elseif(k_t > domhi(3)) then
              write(err_msg, 1101) icv, 'IC_K_T > domhi(3)'
              call flush_err_msg(abort=.true.)
           endif

        endif
     enddo   ! end loop over (icv=1,dim_ic)

1101 format('Error 1101: Initial condition region ',I2,' is ill-',    &
        'defined.',/3x,A,/'Please correct the input deck.')

      call finl_err_msg

   end subroutine check_ic_geometry


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_IC_GAS_PHASE                                       !
!                                                                      !
! Purpose: Verify gas phase input variables in IC region.              !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   subroutine check_ic_gas_phase(ICV)

      use ic,        only: ic_u, ic_v, ic_w

      integer, intent(in) :: ICV

      call init_err_msg("CHECK_IC_GAS_PHASE")

      ! Check that gas phase velocity components are initialized.
      if(is_undefined(ic_u(icv))) then
         write(err_msg, 1000) trim(ivar('ic_u',icv))
         call flush_err_msg(abort=.true.)
      endif

      if(is_undefined(ic_v(icv))) then
         write(err_msg, 1000) trim(ivar('ic_v',icv))
         call flush_err_msg(abort=.true.)
      endif

      if(is_undefined(ic_w(icv))) then
         write(err_msg, 1000) trim(ivar('ic_w',icv))
         call flush_err_msg(abort=.true.)
      endif

1000  format('Error 1000: Required input not specified: ',A,/&
         'Please correct the input deck.')

      call finl_err_msg

   end subroutine check_ic_gas_phase

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_IC_OVERFLOW                                        !
!                                                                      !
! Purpose: Verify that no data was defined for unspecified IC.         !
!                                                                      !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    subroutine check_ic_overflow(ICV)

      use ic,    only: ic_u, ic_v, ic_w

      integer, intent(in) :: icv

      ! Initialize the error manager.
      call init_err_msg("CHECK_IC_OVERFLOW")

      if(is_defined(ic_u(icv))) then
         write(err_msg, 1010) trim(ivar('IC_U',ICV))
         call flush_err_msg(abort=.true.)

      elseif(is_defined(ic_v(icv))) then
         write(err_msg, 1010) trim(ivar('IC_V',ICV))
         call flush_err_msg(abort=.true.)

      elseif(is_defined(ic_w(icv))) then
         write(err_msg, 1010) trim(ivar('IC_W',ICV))
         call flush_err_msg(abort=.true.)

      endif

      call finl_err_msg

1010  format('Error 1010: ',A,' specified in an undefined IC region')

    end subroutine check_ic_overflow

end module check_initial_conditions_module
