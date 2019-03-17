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
                          ng) bind(c,name='set_bc_type')

      use bc,       only: bc_defined, bc_type, bc_plane
      use bc,       only: cyclic_x, cyclic_y, cyclic_z
      use bc,       only: minf_, nsw_, pinf_, pout_
      use bc,       only: undef_cell
      use constant, only: dim_bc

      implicit none

      integer(c_int), intent(in   ) :: ng
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)

      real(rt),       intent(in   ) :: dx, dy, dz
      real(rt)  ,     intent(in   ) :: xlength, ylength, zlength

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

subroutine set_bc_mod(pid, ptype, plo, phi, ploc, ppg, pvel) &
     bind(C,name='set_bc_mod')

  use bc

  implicit none

  integer(c_int), intent(in   ) :: pid, ptype
  real(rt),       intent(in   ) :: plo(3), phi(3), ploc, ppg, pvel(3)

  real(rt), parameter :: offset = 1.0d-15

  select case(pid)

  case(1); bc_plane(pid) = 'E'

     bc_center(pid,1) = ploc + offset
     bc_center(pid,2) = plo(2) + 0.5_rt*(phi(2) - plo(2))
     bc_center(pid,3) = plo(3) + 0.5_rt*(phi(3) - plo(3))

     bc_normal(pid,:) = (/ 1.0d0, 0.0d0, 0.0d0/)

  case(2); bc_plane(pid) = 'W'

     bc_center(pid,1) = ploc - offset
     bc_center(pid,2) = plo(2) + 0.5_rt*(phi(2) - plo(2))
     bc_center(pid,3) = plo(3) + 0.5_rt*(phi(3) - plo(3))

     bc_normal(pid,:) = (/-1.0d0, 0.0d0, 0.0d0/)

  case(3); bc_plane(pid) = 'N'

     bc_center(pid,1) = plo(1) + 0.5_rt*(phi(1) - plo(1))
     bc_center(pid,2) = ploc + offset
     bc_center(pid,3) = plo(3) + 0.5_rt*(phi(3) - plo(3))

     bc_normal(pid,:) = (/ 0.0d0, 1.0d0, 0.0d0/)

  case(4); bc_plane(pid) = 'S'

     bc_center(pid,1) = plo(1) + 0.5_rt*(phi(1) - plo(1))
     bc_center(pid,2) = ploc - offset
     bc_center(pid,3) = plo(3) + 0.5_rt*(phi(3) - plo(3))

     bc_normal(pid,:) = (/ 0.0d0,-1.0d0, 0.0d0/)

  case(5); bc_plane(pid) = 'T'

     bc_center(pid,1) = plo(1) + 0.5_rt*(phi(1) - plo(1))
     bc_center(pid,2) = plo(2) + 0.5_rt*(phi(2) - plo(2))
     bc_center(pid,3) = ploc + offset

     bc_normal(pid,:) = (/ 0.0d0, 0.0d0, 1.0d0/)

  case(6); bc_plane(pid) = 'B'

     bc_center(pid,1) = plo(1) + 0.5_rt*(phi(1) - plo(1))
     bc_center(pid,2) = plo(2) + 0.5_rt*(phi(2) - plo(2))
     bc_center(pid,3) = ploc - offset

     bc_normal(pid,:) = (/ 0.0d0, 0.0d0,-1.0d0/)

  end select


  select case(ptype)

  case(MINF_)

     bc_type(pid) = 'MI'

     bc_p(pid) =   ppg;

     bc_u(pid) = pvel(1);
     bc_v(pid) = pvel(2);
     bc_w(pid) = pvel(3);

     bc_defined(pid) = .true.

  case(PINF_)

     bc_type(pid) = 'PI'

     bc_p(pid) =   ppg;

     bc_defined(pid) = .true.

  case(POUT_)

     bc_type(pid) = 'PO'

     bc_p(pid) =   ppg;

     bc_defined(pid) = .true.

  case(NSW_)

     bc_type(pid) = 'NSW'

     bc_u(pid) = pvel(1)
     bc_v(pid) = pvel(2)
     bc_w(pid) = pvel(3)

     select case(pid)
     case(1,2); bc_u(pid) = 0.0d0;
     case(3,4); bc_v(pid) = 0.0d0;
     case(5,6); bc_w(pid) = 0.0d0;
     end select

     bc_defined(pid) = .true.

  case DEFAULT

     bc_defined(pid) = .false.

  end select


end subroutine set_bc_mod

end module set_bc_type_module
