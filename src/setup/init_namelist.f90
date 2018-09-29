MODULE INIT_NAMELIST_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: INIT_NAMELIST                                           !
!  Purpose: initialize the NAMELIST variables                          !
!                                                                      !
!  Author: P. Nicoletti                               Date: 26-NOV-91  !
!                                                                      !
!  Keyword Documentation Format:                                       !
!                                                                      !
!<keyword category="category name" required="true"/FALSE               !
!                                    legacy=TRUE/FALSE>                !
!  <description></description>                                         !
!  <arg index="" id="" max="" min=""/>                                 !
!  <dependent keyword="" value="DEFINED"/>                             !
!  <conflict keyword="" value="DEFINED"/>                              !
!  <valid value="" note="" alias=""/>                                  !
!  <range min="" max="" />                                             !
!  INCFLO_KEYWORD=INIT_VALUE                                             !
!</keyword>                                                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE INIT_NAMELIST

      use bc
      use ic
      use constant, only: gravity
      use ic, only: ic_p, ic_t
      use ic, only: ic_u, ic_v, ic_w
      use ic, only: ic_x_e, ic_x_w, ic_y_n, ic_y_s, ic_z_b, ic_z_t
      use utilities, only: blank_line, seek_comment
      use utilities, only: make_upper_case, replace_tab

      use param, only: zero, one
      use param, only: undefined, undefined_c

      implicit none

!#####################################################################!
!                      Geometry and Discretization                    !
!#####################################################################!

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Fluid pressure drop across xlength when a cyclic boundary condition
!    with pressure drop is imposed in the x-direction.
!  </description>
      delp_x = zero
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Fluid pressure drop across ylength when a cyclic boundary condition
!    with pressure drop is imposed in the y-direction.
!  </description>
      delp_y = zero
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Fluid pressure drop across zlength when a cyclic boundary condition
!    with pressure drop is imposed in the z-direction.
!  </description>
      delp_z = zero
!</keyword>


!#####################################################################!
!                   Initial Conditions Section                        !
!#####################################################################!


!<keyword category="Initial Condition" required="false">
!  <description>X coordinate of the west face.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_X_W(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>X coordinate of the east face.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_X_E(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Y coordinate of the south face.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_Y_S(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Y coordinate of the north face.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_Y_N(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Z coordinate of the bottom face.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_Z_B(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Z coordinate of the top face.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_Z_T(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Initial gas pressure in the IC region. If this quantity is not
!    specified, incflo will set up a hydrostatic pressure profile,
!    which varies only in the y-direction.
!  </description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_P(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial gas phase temperature in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_T(:) = 293.15d0
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial x-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_U(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial y-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_V(:) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial z-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIM_IC"/>
      IC_W(:) = UNDEFINED
!</keyword>

!#####################################################################!
!                        Boundary Conditions                          !
!#####################################################################!

!<keyword category="Boundary Condition" required="false">
!  <description>X coordinate of the west face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_X_W(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>X coordinate of the east face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_X_E(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Y coordinate of the south face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_Y_S(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Y coordinate of the north face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_Y_N(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z coordinate of the bottom face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_Z_B(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z coordinate of the top face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_Z_T(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z coordinate of the top face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_NORMAL(:,:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z coordinate of the top face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_CENTER(:,:) = UNDEFINED
!</keyword>


!<keyword category="Boundary Condition" required="false">
!  <description>Type of boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
!
!  <valid value='DUMMY'
!    note='The specified boundary condition is ignored. This is
!      useful for turning off some boundary conditions without having
!      to delete them from the file.' />
!
!  <valid value='MASS_INFLOW' alias='MI'
!    note='Mass inflow rates for gas and solids phases are
!      specified at the boundary.'/>
!
!  <valid value='MASS_OUTFLOW' alias='MO'
!    note='The specified values of gas and solids mass outflow
!      rates at the boundary are maintained, approximately. This
!      condition should be used sparingly for minor outflows, when
!      the bulk of the outflow is occurring through other constant
!      pressure outflow boundaries.' />
!
!  <valid value='P_INFLOW' alias='PI'
!    note='Inflow from a boundary at a specified constant
!      pressure. To specify as the west, south, or bottom end of
!      the computational region, add a layer of wall cells to the
!      west, south, or bottom of the PI cells. Users need to specify
!      all scalar quantities and velocity components. The specified
!      values of fluid and solids velocities are only used initially
!      as incflo computes these values at this inlet boundary.' />
!
!  <valid value='FREE_SLIP_WALL' alias='FSW'
!    note='Velocity gradients at the wall vanish./>
!
!  <valid value='NO_SLIP_WALL' alias='NSW'
!    note='All components of the velocity vanish at the wall./>
!
!  <valid value='PAR_SLIP_WALL' alias='PSW'
!    note='Partial slip at the wall implemented as
!      dv/dn + hw (v - vw) = 0, where n is the normal pointing from the
!      fluid into the wall. The coefficients hw and vw should be
!      specified. For free slip set hw = 0. For no slip leave hw
!      undefined (hw=+inf) and set vw = 0. To set hw = +inf, leave it
!      unspecified. />
      BC_TYPE(:) = UNDEFINED_C
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase hw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_HW(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Uw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_UW(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Vw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_VW(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Ww for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_WW(:) = UNDEFINED
!</keyword>


!<keyword category="Boundary Condition" required="false">
!  <description>
!    Gas phase heat transfer coefficient, Hw, in diffusion boundary condition:
!    d(T)/dn + Hw (T - Tw) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_HW_T(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified gas phase wall temperature, Tw, in diffusion boundary condition:
!    d(T)/dn + Hw (T - Tw) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_TW(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified constant gas phase heat flux, C, in diffusion boundary condition:
!    d(T)/dn + Hw (T - Tw) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_C_T(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas pressure at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_P(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase temperature at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_T(:) = 293.15d0
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>X-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_U(:) = ZERO
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Y-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_V(:) = ZERO
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Z-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_W(:) = ZERO
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas volumetric flow rate through the boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_VOLFLOW(:) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas mass flow rate through the boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIM_BC"/>
      BC_MASSFLOW(:) = UNDEFINED
!</keyword>

      end subroutine init_namelist
end module init_namelist_module
