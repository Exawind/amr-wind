module set_bc_flow_module

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int

  use bc, only: bc_type, bc_plane

  use param,  only: dim_m
  use param, only: zero, one, equal, is_defined

  use bc, only: bc_u_g, bc_v_g, bc_w_g
  use bc, only: bc_massflow_g, bc_volflow_g

  use bc, only: bc_u_s, bc_v_s, bc_w_s
  use bc, only: bc_massflow_s, bc_volflow_s

  implicit none
  private 

  public set_bc_flow

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_FLOW                                             !
!                                                                      !
!  Purpose: Check boundary condition specifications                    !
!     - convert physical locations to i, j, k's                        !
!     - convert mass and volumetric flows to velocities (FLOW_TO_VEL)  !
!     - check specification of physical quantities                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine set_bc_flow(xlength, ylength, zlength, dx, dy, dz) &
     bind(C,name ="set_bc_flow")

     use param,    only: dim_bc
     use bc,       only: bc_defined, bc_type

    implicit none

    real(rt)  , intent(in) :: xlength, ylength, zlength
    real(rt)  , intent(in) :: dx, dy, dz

    integer :: bcv

    ! Loop over each defined BC and check the user data.

    do bcv = 1, dim_bc
       if(bc_defined(bcv)) then

          select case (trim(bc_type(bcv)))

             case ('MASS_INFLOW','MI','MASS_OUTFLOW','MO')
                call gas_volflow_to_vel(bcv, xlength, ylength, zlength, dx, dy, dz)

          end select
       endif
    enddo

  end subroutine set_bc_flow

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GAS_VOLFLOW_TO_VELOCITY                                 !
!                                                                      !
!  Purpose: Convert gas phase volumetric rate to a velocity.           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine gas_volflow_to_vel(bcv, xlength, ylength, zlength, &
     dx, dy, dz)

     use bc, only: bc_x_w, bc_y_s, bc_z_b
     use bc, only: bc_x_e, bc_y_n, bc_z_t
     use calc_cell_module, only: calc_cell_bc_flow

     implicit none

     integer,        intent(in) :: bcv
     real(rt)  , intent(in) :: xlength, ylength, zlength
     real(rt)  , intent(in) :: dx, dy, dz

     real(rt) :: sgn, off, vel, area
     integer  :: i_w, i_e, j_s, j_n, k_b, k_t

    select case (trim(bc_type(bcv)))
    case ('MASS_INFLOW', 'MI'); SGN =  ONE; OFF = ZERO
    case ('MASS_OUTFLOW','MO'); SGN = -ONE; OFF = ONE
    end select

    select case (bc_plane(bcv))
    case ('W'); sgn = -sgn
    case ('S'); sgn = -sgn
    case ('B'); sgn = -sgn
    end select

    call calc_cell_bc_flow(&
       xlength, ylength, zlength, dx, dy, dz, &
       bc_x_w(bcv), bc_y_s(bcv), bc_z_b(bcv), &
       bc_x_e(bcv), bc_y_n(bcv), bc_z_t(bcv), &
       i_w, i_e, j_s, j_n, k_b, k_t)

    select case(bc_plane(bcv))
    case('W','E')
       area = dy*dble(j_n-j_s+1)*dz*dble(k_t-k_b+1)
       vel = sgn*bc_volflow_g(bcv)/(area)
       bc_u_g(bcv) = vel
       bc_v_g(bcv) = off * bc_v_g(bcv)
       bc_w_g(bcv) = off * bc_w_g(bcv)
    case('S','N')
       area = dx*dble(i_e-i_w+1)*dz*dble(k_t-k_b+1)
       vel = sgn*bc_volflow_g(bcv)/(area)
       bc_v_g(bcv) = vel
       bc_u_g(bcv) = off * bc_u_g(bcv)
       bc_w_g(bcv) = off * bc_w_g(bcv)
    case('B','T')
       area = dx*dble(i_e-i_w+1)*dy*dble(j_n-j_s+1)
       vel = sgn*bc_volflow_g(bcv)/(area)
       bc_w_g(bcv) = vel
       bc_u_g(bcv) = off * bc_u_g(bcv)
       bc_v_g(bcv) = off * bc_v_g(bcv)
    end select

  end subroutine gas_volflow_to_vel

end module set_bc_flow_module
