!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                       C
!  Subroutine: zero_wall_norm_vel                                       C
!                                                                       C
!  This subroutine does the initial setting of all boundary conditions. C
!                                                                       C
!                                                                       C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
subroutine zero_wall_norm_vel(slo, shi, &
     vel, bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
     bc_klo_type, bc_khi_type, domlo, domhi, ng) &
     bind(C, name="zero_wall_norm_vel")

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int

  use bc, only: nsw_, fsw_, psw_, pinf_, pout_, minf_

  implicit none

  integer(c_int), intent(in   ) :: slo(3),shi(3)
  integer(c_int), intent(in   ) :: domlo(3),domhi(3),ng

  real(rt), intent(inout) :: &
       vel(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

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
  integer    ilo, ihi, jlo, jhi, klo, khi

!--------------------------------------------------------------------//

  nlft = max(0,domlo(1)-slo(1))
  nbot = max(0,domlo(2)-slo(2))
  ndwn = max(0,domlo(3)-slo(3))

  nrgt = max(0,shi(1)-domhi(1))
  ntop = max(0,shi(2)-domhi(2))
  nup  = max(0,shi(3)-domhi(3))

  if (nlft .gt. 0) then
     ilo = domlo(1)
     do i = 1, nlft
        do k=slo(3),shi(3)
           do j=slo(2),shi(2)
              bcv = bc_ilo_type(j,k,2)
              if(bc_ilo_type(j,k,1) == NSW_ .or. &
                 bc_ilo_type(j,k,1) == PSW_ .or. &
                 bc_ilo_type(j,k,1) == FSW_) then

                 vel(ilo-i,j,k,1) = 0.d0

              end if
           end do
        end do
     end do
  endif

  if (nrgt .gt. 0) then
     ihi = domhi(1)
     do i = 0, nrgt
        do k=slo(3),shi(3)
           do j=slo(2),shi(2)
              bcv = bc_ihi_type(j,k,2)
              if(bc_ihi_type(j,k,1) == NSW_ .or. &
                 bc_ihi_type(j,k,1) == PSW_ .or. &
                 bc_ihi_type(j,k,1) == FSW_) then

                 vel(ihi+i,j,k,1) = 0.d0

              end if
           end do
        end do
     end do
  endif

  if (nbot .gt. 0) then
     jlo = domlo(2)
     do j = 1, nbot
        do k=slo(3),shi(3)
           do i=slo(1),shi(1)
              bcv = bc_jlo_type(i,k,2)
              if(bc_jlo_type(i,k,1) == NSW_ .or. &
                 bc_jlo_type(i,k,1) == PSW_ .or. &
                 bc_jlo_type(i,k,1) == FSW_) then

                 vel(i,jlo-j,k,2) = 0.d0

              end if
           end do
        end do
     end do
  endif

  if (ntop .gt. 0) then
     jhi = domhi(2)
     do j = 0, ntop
        do k=slo(3),shi(3)
           do i=slo(1),shi(1)
              bcv = bc_jhi_type(i,k,2)
              if(bc_jhi_type(i,k,1) == NSW_ .or. &
                 bc_jhi_type(i,k,1) == PSW_ .or. &
                 bc_jhi_type(i,k,1) == FSW_) then

                 vel(i,jhi+j,k,2) = 0.d0

              end if
           end do
        end do
     end do
  endif

  if (ndwn .gt. 0) then
     klo = domlo(3)
     do k = 1, ndwn
        do j=slo(2),shi(2)
           do i=slo(1),shi(1)
              bcv = bc_klo_type(i,j,2)
              if(bc_klo_type(i,j,1) == NSW_ .or. &
                 bc_klo_type(i,j,1) == PSW_ .or. &
                 bc_klo_type(i,j,1) == FSW_) then

                 vel(i,j,klo-k,3) = 0.d0

              end if
           end do
        end do
     end do
  endif

  if (nup .gt. 0) then
     khi = domhi(3)
     do k = 0, nup
        do j=slo(2),shi(2)
           do i=slo(1),shi(1)
              bcv = bc_khi_type(i,j,2)
              if(bc_khi_type(i,j,1) == NSW_ .or. &
                 bc_khi_type(i,j,1) == PSW_ .or. &
                 bc_khi_type(i,j,1) == FSW_) then

                 vel(i,j,khi+k,3) = 0.d0

              end if
           end do
        end do
     end do
  endif

end subroutine zero_wall_norm_vel
