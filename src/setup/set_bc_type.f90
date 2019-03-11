module set_bc_type_module

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc_type                                             C
!                                                                      C
!  Author: J. Musser                                  Date: 05-FEB-17  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   subroutine set_bc_type(bc_ilo_type, bc_ihi_type, &
                          bc_jlo_type, bc_jhi_type, &
                          bc_klo_type, bc_khi_type, &
                          domlo, domhi, dx, dy, dz, &
                          xlength, ylength, zlength,&
                          ng) &
               bind(c,name='set_bc_type')

      use bc, only: bc_defined, bc_type, bc_plane

      use bc, only: nsw_, pinf_, pout_, minf_
      use bc, only: undef_cell
      use bc, only: cyclic_x, cyclic_y, cyclic_z

      use param, only: dim_bc
      use param, only: equal
      use calc_cell_module, only: calc_cell_bc_flow
      use calc_cell_module, only: calc_cell_bc_wall

      implicit none

      integer(c_int), intent(in   ) :: ng
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)

      real(rt), intent(in) :: dx, dy, dz
      real(rt)  , intent(in) :: xlength, ylength, zlength

      integer(c_int), intent(inout) :: bc_ilo_type&
         (domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2)
      integer(c_int), intent(inout) :: bc_ihi_type&
         (domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2)
      integer(c_int), intent(inout) :: bc_jlo_type&
         (domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2)
      integer(c_int), intent(inout) :: bc_jhi_type&
         (domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2)
      integer(c_int), intent(inout) :: bc_klo_type&
         (domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)
      integer(c_int), intent(inout) :: bc_khi_type&
         (domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

      ! Local index for boundary condition
      integer :: type, bcv
      integer :: i,j,k

      integer :: i_w, j_s, k_b, i_e, j_n, k_t

      bc_ilo_type(:,:,1) = merge(undef_cell, nsw_, cyclic_x)
      bc_ihi_type(:,:,1) = merge(undef_cell, nsw_, cyclic_x)
      bc_jlo_type(:,:,1) = merge(undef_cell, nsw_, cyclic_y)
      bc_jhi_type(:,:,1) = merge(undef_cell, nsw_, cyclic_y)
      bc_klo_type(:,:,1) = merge(undef_cell, nsw_, cyclic_z)
      bc_khi_type(:,:,1) = merge(undef_cell, nsw_, cyclic_z)

      do bcv = 1, dim_bc
         if (bc_defined(bcv)) then

            select case (trim(bc_type(bcv)))
            case('NO_SLIP_WALL'  ,'NSW'); type = nsw_
            case('P_INFLOW'      ,'PI' ); type = pinf_
            case('P_OUTFLOW'     ,'PO' ); type = pout_
            case('MASS_INFLOW'   ,'MI' ); type = minf_
            case default
               write(6,*) 'unknown bc type'
               stop 7655
            end select

            if (bc_plane(bcv) == 'E') then
               bc_ilo_type(:,:,1) = type
               bc_ilo_type(:,:,2) = bcv

            else if(bc_plane(bcv) == 'W') then
               bc_ihi_type(:,:,1) = type
               bc_ihi_type(:,:,2) = bcv

            else if(bc_plane(bcv) == 'N') then
               bc_jlo_type(:,:,1) = type
               bc_jlo_type(:,:,2) = bcv

            else if(bc_plane(bcv) == 'S') then
               bc_jhi_type(:,:,1) = type
               bc_jhi_type(:,:,2) = bcv

            else if(bc_plane(bcv) == 'T') then
               bc_klo_type(:,:,1) = type
               bc_klo_type(:,:,2) = bcv

            else if(bc_plane(bcv) == 'B') then
               bc_khi_type(:,:,1) = type
               bc_khi_type(:,:,2) = bcv

            endif

         endif
      enddo

   end subroutine set_bc_type

SUBROUTINE SET_BC_MOD(PID, PTYPE, PLO, PHI, PLOC, PPG, PVEL) &
     BIND(C,NAME='set_bc_mod')

  USE BC, ONLY: BC_DEFINED
  USE BC, ONLY: BC_TYPE, BC_PLANE

  USE BC, ONLY: NSW_, PINF_, POUT_, MINF_

  USE BC, ONLY: BC_CENTER
  USE BC, ONLY: BC_NORMAL

  USE BC, ONLY: BC_P
  USE BC, ONLY: BC_U, BC_V, BC_W

  IMPLICIT NONE

  INTEGER(C_INT), INTENT(IN   ) :: PID, PTYPE
  REAL(RT),       INTENT(IN   ) :: PLO(3), PHI(3), PLOC, PPG, PVEL(3)

  REAL(RT), PARAMETER :: OFFSET = 1.0D-15

  SELECT CASE(PID)

  CASE(1); BC_PLANE(PID) = 'E'

     BC_CENTER(PID,1) = PLOC + OFFSET
     BC_CENTER(PID,2) = PLO(2) + 0.5_RT*(PHI(2) - PLO(2))
     BC_CENTER(PID,3) = PLO(3) + 0.5_RT*(PHI(3) - PLO(3))

     BC_NORMAL(PID,:) = (/ 1.0D0, 0.0D0, 0.0D0/)

  CASE(2); BC_PLANE(PID) = 'W'

     BC_CENTER(PID,1) = PLOC - OFFSET
     BC_CENTER(PID,2) = PLO(2) + 0.5_RT*(PHI(2) - PLO(2))
     BC_CENTER(PID,3) = PLO(3) + 0.5_RT*(PHI(3) - PLO(3))

     BC_NORMAL(PID,:) = (/-1.0D0, 0.0D0, 0.0D0/)

  CASE(3); BC_PLANE(PID) = 'N'

     BC_CENTER(PID,1) = PLO(1) + 0.5_RT*(PHI(1) - PLO(1))
     BC_CENTER(PID,2) = PLOC + OFFSET
     BC_CENTER(PID,3) = PLO(3) + 0.5_RT*(PHI(3) - PLO(3))

     BC_NORMAL(PID,:) = (/ 0.0D0, 1.0D0, 0.0D0/)

  CASE(4); BC_PLANE(PID) = 'S'

     BC_CENTER(PID,1) = PLO(1) + 0.5_RT*(PHI(1) - PLO(1))
     BC_CENTER(PID,2) = PLOC - OFFSET
     BC_CENTER(PID,3) = PLO(3) + 0.5_RT*(PHI(3) - PLO(3))

     BC_NORMAL(PID,:) = (/ 0.0D0,-1.0D0, 0.0D0/)

  CASE(5); BC_PLANE(PID) = 'T'

     BC_CENTER(PID,1) = PLO(1) + 0.5_RT*(PHI(1) - PLO(1))
     BC_CENTER(PID,2) = PLO(2) + 0.5_RT*(PHI(2) - PLO(2))
     BC_CENTER(PID,3) = PLOC + OFFSET

     BC_NORMAL(PID,:) = (/ 0.0D0, 0.0D0, 1.0D0/)

  CASE(6); BC_PLANE(PID) = 'B'

     BC_CENTER(PID,1) = PLO(1) + 0.5_RT*(PHI(1) - PLO(1))
     BC_CENTER(PID,2) = PLO(2) + 0.5_RT*(PHI(2) - PLO(2))
     BC_CENTER(PID,3) = PLOC - OFFSET

     BC_NORMAL(PID,:) = (/ 0.0D0, 0.0D0,-1.0D0/)

  END SELECT


  SELECT CASE(PTYPE)

  CASE(MINF_)

     BC_TYPE(PID) = 'MI'

     BC_P(PID) =   PPG;

     BC_U(PID) = PVEL(1);
     BC_V(PID) = PVEL(2);
     BC_W(PID) = PVEL(3);

     BC_DEFINED(PID) = .TRUE.

  CASE(PINF_)

     BC_TYPE(PID) = 'PI'

     BC_P(PID) =   PPG;

     BC_DEFINED(PID) = .TRUE.

  CASE(POUT_)

     BC_TYPE(PID) = 'PO'

     BC_P(PID) =   PPG;

     BC_DEFINED(PID) = .TRUE.

  CASE(NSW_)

     BC_TYPE(PID) = 'NSW'

     BC_U(PID) = PVEL(1)
     BC_V(PID) = PVEL(2)
     BC_W(PID) = PVEL(3)

     SELECT CASE(PID)
     CASE(1,2); BC_U(PID) = 0.0D0;
     CASE(3,4); BC_V(PID) = 0.0D0;
     CASE(5,6); BC_W(PID) = 0.0D0;
     END SELECT

     BC_DEFINED(PID) = .TRUE.

  CASE DEFAULT

     BC_DEFINED(PID) = .FALSE.

  END SELECT


END SUBROUTINE SET_BC_MOD

end module set_bc_type_module
