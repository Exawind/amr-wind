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

      use bc, only: nsw_, fsw_, psw_
      use bc, only: pinf_, pout_
      use bc, only: minf_
      use bc, only: undef_cell, cycl_
      use bc, only: cyclic_x, cyclic_y, cyclic_z

      use bc, only: bc_x_w, bc_y_s, bc_z_b
      use bc, only: bc_x_e, bc_y_n, bc_z_t

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
               case('FREE_SLIP_WALL','FSW'); type = fsw_
               case('NO_SLIP_WALL'  ,'NSW'); type = nsw_
               case('PAR_SLIP_WALL' ,'PSW'); type = psw_
               case('P_INFLOW'      ,'PI' ); type = pinf_
               case('P_OUTFLOW'     ,'PO' ); type = pout_
               case('MASS_INFLOW'   ,'MI' ); type = minf_
               case default
                  write(6,*) 'unknown bc type'
                  stop 7655
            end select

            select case(type)
            case(nsw_, fsw_, psw_)
               call calc_cell_bc_wall(domlo, domhi, &
                  xlength, ylength, zlength, dx, dy, dz, &
                  bc_x_w(bcv), bc_y_s(bcv), bc_z_b(bcv), &
                  bc_x_e(bcv), bc_y_n(bcv), bc_z_t(bcv), &
                  i_w, i_e, j_s, j_n, k_b, k_t)
            case(pinf_, pout_, minf_)
               call calc_cell_bc_flow(&
                  xlength, ylength, zlength, dx, dy, dz, &
                  bc_x_w(bcv), bc_y_s(bcv), bc_z_b(bcv), &
                  bc_x_e(bcv), bc_y_n(bcv), bc_z_t(bcv), &
                  i_w, i_e, j_s, j_n, k_b, k_t)
            end select

            if (i_w == i_e) then
               if(i_w == domlo(1)-1) then
                  bc_ilo_type(j_s:j_n, k_b:k_t,1) = type
                  bc_ilo_type(j_s:j_n, k_b:k_t,2) = bcv

                  bc_plane(bcv) = 'E' ! Fluid is to East of bc!

               else if(i_w == domhi(1)+1) then
                  bc_ihi_type(j_s:j_n, k_b:k_t,1) = type
                  bc_ihi_type(j_s:j_n, k_b:k_t,2) = bcv

                  bc_plane(bcv) = 'W' ! Fluid is to West of bc!

               endif
            endif

            if (j_s == j_n) then
               if(j_s == domlo(2)-1) then
                  bc_jlo_type(i_w:i_e, k_b:k_t,1) = type
                  bc_jlo_type(i_w:i_e, k_b:k_t,2) = bcv

                  bc_plane(bcv) = 'N' ! Fluid is to North of bc!

               else if(j_s == domhi(2)+1) then
                  bc_jhi_type(i_w:i_e, k_b:k_t,1) = type
                  bc_jhi_type(i_w:i_e, k_b:k_t,2) = bcv

                  bc_plane(bcv) = 'S' ! Fluid is to South of bc!

               endif
            endif

            if (k_b == k_t) then
               if(k_b == domlo(3)-1) then
                  bc_klo_type(i_w:i_e, j_s:j_n,1) = type
                  bc_klo_type(i_w:i_e, j_s:j_n,2) = bcv

                  bc_plane(bcv) = 'T' ! Fluid is to Top of bc!

               elseif(k_b == domhi(3)+1) then
                  bc_khi_type(i_w:i_e, j_s:j_n,1) = type
                  bc_khi_type(i_w:i_e, j_s:j_n,2) = bcv

                  bc_plane(bcv) = 'B' ! Fluid is to Bottom of bc!

               endif
            endif

         endif
      enddo


! Edge cases
! --------------------------------------------------------------------------------------------

      do j=1,ng
         bc_ilo_type(domlo(2)-j,domlo(3)-1:domhi(3)+1,:) = bc_ilo_type(domlo(2),domlo(3)-1:domhi(3)+1,:)
         bc_ihi_type(domlo(2)-j,domlo(3)-1:domhi(3)+1,:) = bc_ihi_type(domlo(2),domlo(3)-1:domhi(3)+1,:)
      enddo
      do j=1,ng
         bc_ilo_type(domhi(2)+j,domlo(3)-1:domhi(3)+1,:) = bc_ilo_type(domhi(2),domlo(3)-1:domhi(3)+1,:)
         bc_ihi_type(domhi(2)+j,domlo(3)-1:domhi(3)+1,:) = bc_ihi_type(domhi(2),domlo(3)-1:domhi(3)+1,:)
      enddo
      do k=1,ng
         bc_ilo_type(domlo(2)-1:domhi(2)+1,domlo(3)-k,:) = bc_ilo_type(domlo(2)-1:domhi(2)+1,domlo(3),:)
         bc_ihi_type(domlo(2)-1:domhi(2)+1,domlo(3)-k,:) = bc_ihi_type(domlo(2)-1:domhi(2)+1,domlo(3),:)
      enddo
      do k=1,ng
         bc_ilo_type(domlo(2)-1:domhi(2)+1,domhi(3)+k,:) = bc_ilo_type(domlo(2)-1:domhi(2)+1,domhi(3),:)
         bc_ihi_type(domlo(2)-1:domhi(2)+1,domhi(3)+k,:) = bc_ihi_type(domlo(2)-1:domhi(2)+1,domhi(3),:)
      enddo


      do i=1,ng
         bc_jlo_type(domlo(1)-i,domlo(3)-1:domhi(3)+1,:) = bc_jlo_type(domlo(1),domlo(3)-1:domhi(3)+1,:)
         bc_jhi_type(domlo(1)-i,domlo(3)-1:domhi(3)+1,:) = bc_jhi_type(domlo(1),domlo(3)-1:domhi(3)+1,:)
      enddo
      do i=1,ng
         bc_jlo_type(domhi(1)+i,domlo(3)-1:domhi(3)+1,:) = bc_jlo_type(domhi(1),domlo(3)-1:domhi(3)+1,:)
         bc_jhi_type(domhi(1)+i,domlo(3)-1:domhi(3)+1,:) = bc_jhi_type(domhi(1),domlo(3)-1:domhi(3)+1,:)
      enddo
      do k=1,ng
         bc_jlo_type(domlo(1)-1:domhi(1)+1,domlo(3)-k,:) = bc_jlo_type(domlo(1)-1:domhi(1)+1,domlo(3),:)
         bc_jhi_type(domlo(1)-1:domhi(1)+1,domlo(3)-k,:) = bc_jhi_type(domlo(1)-1:domhi(1)+1,domlo(3),:)
      enddo
      do k=1,ng
         bc_jlo_type(domlo(1)-1:domhi(1)+1,domhi(3)+k,:) = bc_jlo_type(domlo(1)-1:domhi(1)+1,domhi(3),:)
         bc_jhi_type(domlo(1)-1:domhi(1)+1,domhi(3)+k,:) = bc_jhi_type(domlo(1)-1:domhi(1)+1,domhi(3),:)
      enddo


      do i=1,ng
         bc_klo_type(domlo(1)-i,domlo(2)-1:domhi(2)+1,:) = bc_klo_type(domlo(1),domlo(2)-1:domhi(2)+1,:)
         bc_khi_type(domlo(1)-i,domlo(2)-1:domhi(2)+1,:) = bc_khi_type(domlo(1),domlo(2)-1:domhi(2)+1,:)
      enddo
      do i=1,ng
         bc_klo_type(domhi(1)+i,domlo(2)-1:domhi(2)+1,:) = bc_klo_type(domhi(1),domlo(2)-1:domhi(2)+1,:)
         bc_khi_type(domhi(1)+i,domlo(2)-1:domhi(2)+1,:) = bc_khi_type(domhi(1),domlo(2)-1:domhi(2)+1,:)
      enddo
      do j=1,ng
         bc_klo_type(domlo(1)-1:domhi(1)+1,domlo(2)-j,:) = bc_klo_type(domlo(1)-1:domhi(1)+1,domlo(2),:)
         bc_khi_type(domlo(1)-1:domhi(1)+1,domlo(2)-j,:) = bc_khi_type(domlo(1)-1:domhi(1)+1,domlo(2),:)
      enddo
      do j=1,ng
         bc_klo_type(domlo(1)-1:domhi(1)+1,domhi(2)+j,:) = bc_klo_type(domlo(1)-1:domhi(1)+1,domhi(2),:)
         bc_khi_type(domlo(1)-1:domhi(1)+1,domhi(2)+j,:) = bc_khi_type(domlo(1)-1:domhi(1)+1,domhi(2),:)
      enddo

   end subroutine set_bc_type
 end module set_bc_type_module
