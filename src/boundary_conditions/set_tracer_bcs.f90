!
!  This subroutine sets the BCs for the velocity components only.
!
subroutine set_tracer_bcs(time, &
                          tracer, tlo, thi, &
                          bct_ilo, bct_ihi, &
                          bct_jlo, bct_jhi, &
                          bct_klo, bct_khi, &
                          domlo, domhi, &
                          ng, extrap_dir_bcs, probtype) bind(C)

   use amrex_fort_module,  only: ar => amrex_real
   use iso_c_binding,      only: c_int

   use bc
   use constant,           only: zero, half, one, two

   implicit none

   ! Time (necessary if we have time-dependent boundary conditions)
   real(ar),       intent(in   ) :: time

   ! Array bounds
   integer(c_int), intent(in   ) :: tlo(3), thi(3)

   ! This flag, if true, extrapolates Dirichlet bc's to the ghost cells rather than
   !    the faces.
   integer(c_int), intent(in   ) :: extrap_dir_bcs

   ! Grid bounds
   integer(c_int), intent(in   ) :: domlo(3), domhi(3)
   integer(c_int), intent(in   ) :: ng, probtype

   ! BCs type
   integer(c_int), intent(in   )  ::                                 &
        & bct_ilo(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
        & bct_ihi(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
        & bct_jlo(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
        & bct_jhi(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
        & bct_klo(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
        & bct_khi(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

   ! Arrays
   real(ar),      intent(inout) ::  &
      tracer(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3))

   ! Local variables
   integer :: bcv, i, j, k
   integer :: nlft, nrgt, nbot, ntop, nup, ndwn

   ! Used for probtype = {31,32,33} (channel_cylinder with Poiseuille plane inflow BC)
   real    :: x, y, z

   nlft = max(0,domlo(1)-tlo(1))
   nbot = max(0,domlo(2)-tlo(2))
   ndwn = max(0,domlo(3)-tlo(3))

   nrgt = max(0,thi(1)-domhi(1))
   ntop = max(0,thi(2)-domhi(2))
   nup  = max(0,thi(3)-domhi(3))

! *******************************************************************
! First do just the pinf, pout and minf
! *******************************************************************

   if (nlft .gt. 0) then
      do k = tlo(3), thi(3)
         do j = tlo(2), thi(2)

            bcv = bct_ilo(j,k,2)

            select case (bct_ilo(j,k,1))

            case ( pinf_, pout_)

               tracer(tlo(1):domlo(1)-1,j,k) =  tracer(domlo(1),j,k)

            case ( minf_)

               tracer(tlo(1):domlo(1)-1,j,k) = bc_t(bcv)

            end select

         end do
      end do
   endif

   if (nrgt .gt. 0) then

      do k = tlo(3),thi(3)
         do j = tlo(2),thi(2)

            bcv = bct_ihi(j,k,2)

            select case ( bct_ihi(j,k,1) )

            case ( pinf_, pout_ )

               tracer(domhi(1)+1:thi(1),j,k) = tracer(domhi(1),j,k)

            case ( minf_ )

               tracer(domhi(1)+1:thi(1),j,k) = bc_t(bcv)

            end select

         end do
      end do
   endif

   if (nbot .gt. 0) then

      do k = tlo(3), thi(3)
         do i = tlo(1), thi(1)

            bcv = bct_jlo(i,k,2)

            select case ( bct_jlo(i,k,1) )

            case ( pinf_, pout_)

               tracer(i,tlo(2):domlo(2)-1,k) = tracer(i,domlo(2),k)

            case ( minf_ )

               tracer(i,tlo(2):domlo(2)-1,k) = bc_t(bcv)

            end select

         end do
      end do
   endif

   if (ntop .gt. 0) then

      do k = tlo(3), thi(3)
         do i = tlo(1), thi(1)

            bcv = bct_jhi(i,k,2)

            select case ( bct_jhi(i,k,1) )

            case ( pinf_, pout_ )

               tracer(i,domhi(2)+1:thi(2),k) = tracer(i,domhi(2),k)

            case ( minf_)

               tracer(i,domhi(2)+1:thi(2),k) = bc_t(bcv)

            end select

         end do
      end do
   endif

   if (ndwn .gt. 0) then

      do j = tlo(2), thi(2)
         do i = tlo(1), thi(1)

            bcv = bct_klo(i,j,2)

            select case (bct_klo(i,j,1))

            case ( pinf_, pout_ )

               tracer(i,j,tlo(3):domlo(3)-1) = tracer(i,j,domlo(3))

            case ( minf_ )

               tracer(i,j,tlo(3):domlo(3)-1) = bc_t(bcv)

            end select

         end do
      end do
   endif

   if (nup .gt. 0) then

      do j = tlo(2), thi(2)
         do i = tlo(1), thi(1)

            bcv = bct_khi(i,j,2)

            select case ( bct_khi(i,j,1) )

            case ( pinf_, pout_ )

               tracer(i,j,domhi(3)+1:thi(3)) = tracer(i,j,domhi(3))

            case ( minf_ )

               tracer(i,j,domhi(3)+1:thi(3)) = bc_t(bcv)

            end select

         end do
      end do
   endif

! *******************************************************************
! Do this section next to make sure nsw over-rides any previous minf
! *******************************************************************

   if (nlft .gt. 0) then
      do k = tlo(3), thi(3)
         do j = tlo(2), thi(2)

            bcv = bct_ilo(j,k,2)

            select case (bct_ilo(j,k,1))

            case ( nsw_)

               tracer(tlo(1):domlo(1)-1,j,k) = bc_t(bcv)

            end select

         end do
      end do
   endif

   if (nrgt .gt. 0) then

      do k = tlo(3),thi(3)
         do j = tlo(2),thi(2)

            bcv = bct_ihi(j,k,2)

            select case ( bct_ihi(j,k,1) )

            case ( nsw_ )

               tracer(domhi(1)+1:thi(1),j,k) = bc_t(bcv)

            end select

         end do
      end do
   endif

   if (nbot .gt. 0) then

      do k = tlo(3), thi(3)
         do i = tlo(1), thi(1)

            bcv = bct_jlo(i,k,2)

            select case ( bct_jlo(i,k,1) )

            case ( nsw_ )

               tracer(i,tlo(2):domlo(2)-1,k) = bc_t(bcv)

            end select

         end do
      end do
   endif

   if (ntop .gt. 0) then

      do k = tlo(3), thi(3)
         do i = tlo(1), thi(1)

            bcv = bct_jhi(i,k,2)

            select case ( bct_jhi(i,k,1) )

            case ( nsw_)

               tracer(i,domhi(2)+1:thi(2),k) = bc_t(bcv)

            end select

         end do
      end do
   endif

   if (ndwn .gt. 0) then

      do j = tlo(2), thi(2)
         do i = tlo(1), thi(1)

            bcv = bct_klo(i,j,2)

            select case (bct_klo(i,j,1))

            case ( nsw_ )

               tracer(i,j,tlo(3):domlo(3)-1) = bc_t(bcv)

            end select

         end do
      end do
   endif

   if (nup .gt. 0) then

      do j = tlo(2), thi(2)
         do i = tlo(1), thi(1)

            bcv = bct_khi(i,j,2)

            select case ( bct_khi(i,j,1) )

            case ( nsw_ )

               tracer(i,j,domhi(3)+1:thi(3)) = bc_t(bcv)

            end select

         end do
      end do
   endif

end subroutine set_tracer_bcs

