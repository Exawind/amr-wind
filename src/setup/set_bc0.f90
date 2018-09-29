!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc0                                                 C
!  This subroutine does the initial setting of all boundary conditions C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   subroutine set_bc0(slo, shi, &
                      ro, mu, lambda, &
                      bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
                      bc_klo_type, bc_khi_type, domlo, domhi, ng, nodal_pressure &
                      ) bind(C, name="set_bc0")

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      use bc, only: bc_t
      use bc, only: pinf_, pout_, minf_
      use constant, only: ro_0, mu_0

      use param , only: is_undefined

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)
      integer(c_int), intent(in   ) :: ng, nodal_pressure

      real(rt), intent(inout) :: ro&
           (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(inout) :: mu&
           (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(inout) :: lambda&
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

      real(rt) :: bc_ro, bc_mu, bc_lambda
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

                  bc_ro = ro_0

                  bc_mu     = mu_0
                  bc_lambda = -(2.0d0/3.0d0) * mu_0

                      ro(slo(1):domlo(1)-1,j,k) = bc_ro
                      mu(slo(1):domlo(1)-1,j,k) = bc_mu
                  lambda(slo(1):domlo(1)-1,j,k) = bc_lambda

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

                   bc_ro = ro_0

                   bc_mu     = mu_0
                   bc_lambda = -(2.0d0/3.0d0) * mu_0

                        ro(domhi(1)+1:shi(1),j,k) = bc_ro
                        mu(domhi(1)+1:shi(1),j,k) = bc_mu
                    lambda(domhi(1)+1:shi(1),j,k) = bc_lambda

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

                   bc_ro = ro_0

                   bc_mu     = mu_0
                   bc_lambda = -(2.0d0/3.0d0) * mu_0

                      ro(i,slo(2):domlo(2)-1,k) = bc_ro
                      mu(i,slo(2):domlo(2)-1,k) = bc_mu
                  lambda(i,slo(2):domlo(2)-1,k) = bc_lambda

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

                   bc_ro = ro_0

                   bc_mu     = mu_0
                   bc_lambda = -(2.0d0/3.0d0) * mu_0

                      ro(i,domhi(2)+1:shi(2),k) = bc_ro
                      mu(i,domhi(2)+1:shi(2),k) = bc_mu
                  lambda(i,domhi(2)+1:shi(2),k) = bc_lambda

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

                   bc_ro = ro_0

                   bc_mu     = mu_0
                   bc_lambda = -(2.0d0/3.0d0) * mu_0

                       ro(i,j,slo(3):domlo(3)-1) = bc_ro
                       mu(i,j,slo(3):domlo(3)-1) = bc_mu
                   lambda(i,j,slo(3):domlo(3)-1) = bc_lambda

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

                   bc_ro = ro_0

                   bc_mu     = mu_0
                   bc_lambda = -(2.0d0/3.0d0) * mu_0

                       ro(i,j,domhi(3)+1:shi(3)) = bc_ro
                       mu(i,j,domhi(3)+1:shi(3)) = bc_mu
                   lambda(i,j,domhi(3)+1:shi(3)) = bc_lambda

               end if

            end do
         end do
      endif

   end subroutine set_bc0
