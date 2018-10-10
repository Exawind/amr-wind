!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc0                                                 C
!  This subroutine sets initial values of scalar boundary conditions   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   subroutine set_bc0(slo, shi, &
                      ro, eta, &
                      bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
                      bc_klo_type, bc_khi_type, domlo, domhi, ng &
                      ) bind(C, name="set_bc0")

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      use bc, only: bc_t
      use bc, only: pinf_, pout_, minf_
      use constant, only: ro_0, mu

      use param , only: is_undefined

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)
      integer(c_int), intent(in   ) :: ng

      real(rt), intent(inout) :: ro&
                                 (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(inout) :: eta&
                                 (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer(c_int), intent(in   ) :: &
         bc_ilo_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
         bc_ihi_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
         bc_jlo_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
         bc_jhi_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
         bc_klo_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
         bc_khi_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

! Local variables
!--------------------------------------------------------------------//
! local index for boundary condition
      integer :: bcv, i,j,k

      integer    nlft, nrgt, nbot, ntop, nup, ndwn
!--------------------------------------------------------------------//

      nlft = max(0,domlo(1)-slo(1))
      nbot = max(0,domlo(2)-slo(2))
      ndwn = max(0,domlo(3)-slo(3))

      nrgt = max(0,shi(1)-domhi(1))
      ntop = max(0,shi(2)-domhi(2))
      nup  = max(0,shi(3)-domhi(3))

      if (nlft .gt. 0) then
         do k=slo(3),shi(3)
            do j=slo(2),shi(2)

               bcv = bc_ilo_type(j,k,2)

               if (bc_ilo_type(j,k,1) == PINF_ .or. &
                   bc_ilo_type(j,k,1) == POUT_ .or. &
                   bc_ilo_type(j,k,1) == MINF_) then

                  ro(slo(1):domlo(1)-1,j,k) = ro_0
                  eta(slo(1):domlo(1)-1,j,k) = mu

               end if

            end do
         end do
      endif

      if (nrgt .gt. 0) then
         do k=slo(3),shi(3)
            do j=slo(2),shi(2)

               bcv = bc_ihi_type(j,k,2)

               if (bc_ihi_type(j,k,1) == PINF_ .or. &
                   bc_ihi_type(j,k,1) == POUT_ .or. &
                   bc_ihi_type(j,k,1) == MINF_) then

                  ro(domhi(1)+1:shi(1),j,k) = ro_0
                  eta(domhi(1)+1:shi(1),j,k) = mu

               end if

            end do
         end do
      endif

      if (nbot .gt. 0) then
         do k=slo(3),shi(3)
            do i=slo(1),shi(1)

               bcv = bc_jlo_type(i,k,2)

               if (bc_jlo_type(i,k,1) == PINF_ .or. &
                   bc_jlo_type(i,k,1) == POUT_ .or. &
                   bc_jlo_type(i,k,1) == MINF_) then

                  ro(i,slo(2):domlo(2)-1,k) = ro_0
                  eta(i,slo(2):domlo(2)-1,k) = mu

               end if

            end do
         end do
      endif

      if (ntop .gt. 0) then
         do k = slo(3),shi(3)
            do i = slo(1),shi(1)

               bcv = bc_jhi_type(i,k,2)

               if (bc_jhi_type(i,k,1) == PINF_ .or. &
                   bc_jhi_type(i,k,1) == POUT_ .or. &
                   bc_jhi_type(i,k,1) == MINF_) then

                  ro(i,domhi(2)+1:shi(2),k) = ro_0
                  eta(i,domhi(2)+1:shi(2),k) = mu

               end if

            end do
         end do
      endif

      if (ndwn .gt. 0) then
         do j=slo(2),shi(2)
            do i=slo(1),shi(1)

               bcv = bc_klo_type(i,j,2)

               if (bc_klo_type(i,j,1) == PINF_ .or. &
                   bc_klo_type(i,j,1) == POUT_ .or. &
                   bc_klo_type(i,j,1) == MINF_) then

                  ro(i,j,slo(3):domlo(3)-1) = ro_0
                  eta(i,j,slo(3):domlo(3)-1) = mu

               end if

            end do
         end do
      endif

      if (nup .gt. 0) then
         do j=slo(2),shi(2)
            do i=slo(1),shi(1)

               bcv = bc_khi_type(i,j,2)

               if (bc_khi_type(i,j,1) == PINF_ .or. &
                   bc_khi_type(i,j,1) == POUT_ .or. &
                   bc_khi_type(i,j,1) == MINF_) then

                  ro(i,j,domhi(3)+1:shi(3)) = ro_0
                  eta(i,j,domhi(3)+1:shi(3)) = mu

               end if

            end do
         end do
      endif

   end subroutine set_bc0
