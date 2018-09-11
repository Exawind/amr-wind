module check_bc_inflow_module

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int
  use param,         only: undefined, one, zero, is_undefined, is_defined
  use error_manager,  only: finl_err_msg, flush_err_msg, init_err_msg,    &
                         &  ivar,  ival, err_msg

  implicit none
  private

  public check_bc_p_inflow

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: check_bc_p_inflow                                        !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided detailed error message on bc                       !
!                                                                      !
! Comments:                                                            !
!     Unlike the MI boundary, for the PI boundary the velocities at    !
!     the inflow face are calculated by solving the momentum eqns      !
!     and are not fixed. In this way, the PI is similar to the PO      !
!     except that the flow is into the domain and hence all other      !
!     scalars (e.g., mass fractions, void fraction, temperature,       !
!     etc.,) at the inflow cells need to be specified. To satisfy      !
!     the error routines at the start of the simulation, both the      !
!     tangential and normal components at the inflow also need to      !
!     be specified. The velocities values essentially serve as IC.     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  subroutine check_bc_p_inflow(m_tot, skip, bcv)

    use param    , only: dim_m
    use bc,        only: bc_p_g

    integer, intent(in) :: bcv, m_tot
    logical, intent(in) :: skip(dim_m)

    call init_err_msg("CHECK_BC_P_INFLOW")

    if (is_undefined(bc_p_g(bcv))) then
       write(err_msg,1000) 'BC_P_g', bcv
       call flush_err_msg(abort=.true.)
    endif

1000  format('Error 1000: Required input not specified: ',A,/&
             'Please correct the input deck.')

    call finl_err_msg

  end subroutine check_bc_p_inflow

end module check_bc_inflow_module
