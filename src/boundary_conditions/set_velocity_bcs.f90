!
!  This subroutine sets the BCs for the velocity components only.
!
subroutine set_velocity_bcs(time, &
                            vel, ulo, uhi, &
                            bct_ilo, bct_ihi, &
                            bct_jlo, bct_jhi, &
                            bct_klo, bct_khi, &
                            domlo, domhi, &
                            ng, extrap_dir_bcs ) bind(C)

   use amrex_fort_module,  only: ar => amrex_real
   use iso_c_binding    ,  only: c_int
   use param            ,  only: zero,two
   use bc

   implicit none

   ! Time (necessary if we have time-dependent boundary conditions)
   real(ar),       intent(in   ) :: time

   ! Array bounds
   integer(c_int), intent(in   ) :: ulo(3), uhi(3)

   ! This flag, if true, extrapolates Dirichlet bc's to the ghost cells rather than
   !    the faces.
   integer(c_int), intent(in   ) :: extrap_dir_bcs

   ! Grid bounds
   integer(c_int), intent(in   ) :: domlo(3), domhi(3)
   integer(c_int), intent(in   ) :: ng

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
      vel(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)

   ! Local variables
   integer :: bcv, i, j, k
   integer :: nlft, nrgt, nbot, ntop, nup, ndwn

   nlft = max(0,domlo(1)-ulo(1))
   nbot = max(0,domlo(2)-ulo(2))
   ndwn = max(0,domlo(3)-ulo(3))

   nrgt = max(0,uhi(1)-domhi(1))
   ntop = max(0,uhi(2)-domhi(2))
   nup  = max(0,uhi(3)-domhi(3))

   if (nlft .gt. 0) then
      do k = ulo(3), uhi(3)
         do j = ulo(2), uhi(2)

            bcv = bct_ilo(j,k,2)

            select case (bct_ilo(j,k,1))

            case ( pinf_, pout_)

               vel(ulo(1):domlo(1)-1,j,k,1) =  vel(domlo(1),j,k,1)
               vel(ulo(1):domlo(1)-1,j,k,2) =  vel(domlo(1),j,k,2)
               vel(ulo(1):domlo(1)-1,j,k,3) =  vel(domlo(1),j,k,3)

            case ( minf_)

               vel(ulo(1):domlo(1)-1,j,k,1) =  bc_u(bcv)
               vel(ulo(1):domlo(1)-1,j,k,2) =  zero
               vel(ulo(1):domlo(1)-1,j,k,3) =  zero

            case ( nsw_)

               vel(ulo(1):domlo(1)-1,j,k,1) = zero
               vel(ulo(1):domlo(1)-1,j,k,2) = zero
               vel(ulo(1):domlo(1)-1,j,k,3) = zero

            case ( fsw_)

               vel(ulo(1):domlo(1)-1,j,k,1) = zero
               vel(ulo(1):domlo(1)-1,j,k,2) = vel(domlo(1),j,k,2)
               vel(ulo(1):domlo(1)-1,j,k,3) = vel(domlo(1),j,k,3)

            end select

            if (extrap_dir_bcs .gt. 0) then
               select case (bct_ilo(j,k,1))
               case ( minf_, nsw_, fsw_, psw_)
                  vel(domlo(1)-1,j,k,1:3) = two*vel(domlo(1)-1,j,k,1:3) - vel(domlo(1),j,k,1:3)
               end select
            end if

         end do
      end do
   endif

   if (nrgt .gt. 0) then

      do k = ulo(3),uhi(3)
         do j = ulo(2),uhi(2)

            bcv = bct_ihi(j,k,2)

            select case ( bct_ihi(j,k,1) )

            case ( pinf_, pout_ )

               vel(domhi(1)+1:uhi(1),j,k,1) =    vel(domhi(1),j,k,1)
               vel(domhi(1)+1:uhi(1),j,k,2) =    vel(domhi(1),j,k,2)
               vel(domhi(1)+1:uhi(1),j,k,3) =    vel(domhi(1),j,k,3)

            case ( minf_ )

               vel(domhi(1)+1:uhi(1),j,k,1) = bc_u(bcv)
               vel(domhi(1)+1:uhi(1),j,k,2) = zero
               vel(domhi(1)+1:uhi(1),j,k,3) = zero

            case ( nsw_ )

               vel(domhi(1)+1:uhi(1),j,k,1) = zero
               vel(domhi(1)+1:uhi(1),j,k,2) = zero
               vel(domhi(1)+1:uhi(1),j,k,3) = zero

            case ( fsw_ )

               vel(domhi(1)+1:uhi(1),j,k,1) = zero
               vel(domhi(1)+1:uhi(1),j,k,2) = vel(domhi(1),j,k,2)
               vel(domhi(1)+1:uhi(1),j,k,3) = vel(domhi(1),j,k,3)

            end select

            if (extrap_dir_bcs .gt. 0) then
               select case (bct_ihi(j,k,1))
               case ( minf_, nsw_, fsw_, psw_)
                  vel(domhi(1)+1,j,k,1:3) = two*vel(domhi(1)+1,j,k,1:3) - vel(domhi(1),j,k,1:3)
               end select
            end if

         end do
      end do
   endif

   if (nbot .gt. 0) then

      do k = ulo(3), uhi(3)
         do i = ulo(1), uhi(1)

            bcv = bct_jlo(i,k,2)

            select case ( bct_jlo(i,k,1) )

            case ( pinf_, pout_)

               vel(i,ulo(2):domlo(2)-1,k,1) =      vel(i,domlo(2),k,1)
               vel(i,ulo(2):domlo(2)-1,k,2) =      vel(i,domlo(2),k,2)
               vel(i,ulo(2):domlo(2)-1,k,3) =      vel(i,domlo(2),k,3)

            case ( minf_ )

               vel(i,ulo(2):domlo(2)-1,k,1) = zero
               vel(i,ulo(2):domlo(2)-1,k,2) = bc_v(bcv)
               vel(i,ulo(2):domlo(2)-1,k,3) = zero

            case ( nsw_ )

               vel(i,ulo(2):domlo(2)-1,k,1) = zero
               vel(i,ulo(2):domlo(2)-1,k,2) = zero
               vel(i,ulo(2):domlo(2)-1,k,3) = zero

            case ( fsw_)

               vel(i,ulo(2):domlo(2)-1,k,1) = vel(i,domlo(2),k,1)
               vel(i,ulo(2):domlo(2)-1,k,2) = zero
               vel(i,ulo(2):domlo(2)-1,k,3) = vel(i,domlo(2),k,3)

            end select

            if (extrap_dir_bcs .gt. 0) then
               select case (bct_jlo(i,k,1))
               case ( minf_, nsw_, fsw_, psw_)
                  vel(i,domlo(2)-1,k,1:3) = two*vel(i,domlo(2)-1,k,1:3) - vel(i,domlo(2),k,1:3)
               end select
            end if

         end do
      end do
   endif

   if (ntop .gt. 0) then

      do k = ulo(3), uhi(3)
         do i = ulo(1), uhi(1)

            bcv = bct_jhi(i,k,2)

            select case ( bct_jhi(i,k,1) )

            case ( pinf_, pout_ )

               vel(i,domhi(2)+1:uhi(2),k,1) =      vel(i,domhi(2),k,1)
               vel(i,domhi(2)+1:uhi(2),k,2) =      vel(i,domhi(2),k,2)
               vel(i,domhi(2)+1:uhi(2),k,3) =      vel(i,domhi(2),k,3)

            case ( minf_)

               vel(i,domhi(2)+1:uhi(2),k,1) = zero
               vel(i,domhi(2)+1:uhi(2),k,2) = bc_v(bcv)
               vel(i,domhi(2)+1:uhi(2),k,3) = zero

            case ( nsw_)

               vel(i,domhi(2)+1:uhi(2),k,1) = zero
               vel(i,domhi(2)+1:uhi(2),k,2) = zero
               vel(i,domhi(2)+1:uhi(2),k,3) = zero

            case ( fsw_)

               vel(i,domhi(2)+1:uhi(2),k,1) = vel(i,domhi(2),k,1)
               vel(i,domhi(2)+1:uhi(2),k,2) = zero
               vel(i,domhi(2)+1:uhi(2),k,3) = vel(i,domhi(2),k,3)

            end select

            if (extrap_dir_bcs .gt. 0) then
               select case (bct_jhi(i,k,1))
               case ( minf_, nsw_, fsw_, psw_)
                  vel(i,domhi(2)+1,k,1:3) = two*vel(i,domhi(2)+1,k,1:3) - vel(i,domhi(2),k,1:3)
               end select
            end if

         end do
      end do
   endif

   if (ndwn .gt. 0) then

      do j = ulo(2), uhi(2)
         do i = ulo(1), uhi(1)

            bcv = bct_klo(i,j,2)

            select case (bct_klo(i,j,1))

            case ( pinf_, pout_ )

               vel(i,j,ulo(3):domlo(3)-1,1) = vel(i,j,domlo(3),1)
               vel(i,j,ulo(3):domlo(3)-1,2) = vel(i,j,domlo(3),2)
               vel(i,j,ulo(3):domlo(3)-1,3) = vel(i,j,domlo(3),3)

            case ( minf_ )

               vel(i,j,ulo(3):domlo(3)-1,1) = zero
               vel(i,j,ulo(3):domlo(3)-1,2) = zero
               vel(i,j,ulo(3):domlo(3)-1,3) = bc_w(bcv)

            case ( nsw_ )

               vel(i,j,ulo(3):domlo(3)-1,1) = zero
               vel(i,j,ulo(3):domlo(3)-1,2) = zero
               vel(i,j,ulo(3):domlo(3)-1,3) = zero

            case ( fsw_ )

               vel(i,j,ulo(3):domlo(3)-1,1) = vel(i,j,domlo(3),1)
               vel(i,j,ulo(3):domlo(3)-1,2) = vel(i,j,domlo(3),2)
               vel(i,j,ulo(3):domlo(3)-1,3) = zero

            end select

            if (extrap_dir_bcs .gt. 0) then
               select case (bct_klo(i,j,1))
               case ( minf_, nsw_, fsw_, psw_)
                  vel(i,j,domlo(3)-1,1:3) = two*vel(i,j,domlo(3)-1,1:3) - vel(i,j,domlo(3),1:3)
               end select
            end if

         end do
      end do
   endif

   if (nup .gt. 0) then

      do j = ulo(2), uhi(2)
         do i = ulo(1), uhi(1)

            bcv = bct_khi(i,j,2)

            select case ( bct_khi(i,j,1) )

            case ( pinf_, pout_ )

               vel(i,j,domhi(3)+1:uhi(3),1) = vel(i,j,domhi(3),1)
               vel(i,j,domhi(3)+1:uhi(3),2) = vel(i,j,domhi(3),2)
               vel(i,j,domhi(3)+1:uhi(3),3) = vel(i,j,domhi(3),3)

            case ( minf_ )

               vel(i,j,domhi(3)+1:uhi(3),1) = zero
               vel(i,j,domhi(3)+1:uhi(3),2) = zero
               vel(i,j,domhi(3)+1:uhi(3),3) = bc_w(bcv)

            case ( nsw_ )

               vel(i,j,domhi(3)+1:uhi(3),1) = zero
               vel(i,j,domhi(3)+1:uhi(3),2) = zero
               vel(i,j,domhi(3)+1:uhi(3),3) = zero

            case ( fsw_ )

               vel(i,j,domhi(3)+1:uhi(3),1) = vel(i,j,domhi(3),1)
               vel(i,j,domhi(3)+1:uhi(3),2) = vel(i,j,domhi(3),2)
               vel(i,j,domhi(3)+1:uhi(3),3) = zero

            end select

            if (extrap_dir_bcs .gt. 0) then
               select case (bct_khi(i,j,1))
               case ( minf_, nsw_, fsw_, psw_)
                  vel(i,j,domhi(3)+1,1:3) = two*vel(i,j,domhi(3)+1,1:3) - vel(i,j,domhi(3),1:3)
               end select
            end if

         end do
      end do
   endif

   ! *********************************************************************************
   ! We have to do the psw_ bc's last because otherwise non-zero moving wall values
   ! can get over-written
   ! *********************************************************************************

   ! HACK HACK HACK -- we still need to fix the cases below for extrap_dir_bcs = 1

   if (nlft .gt. 0) then
      do k=ulo(3),uhi(3)
         do j=ulo(2),uhi(2)
            bcv = bct_ilo(j,k,2)
            if (bct_ilo(j,k,1) == psw_) then
               vel(ulo(1):domlo(1)-1,j,k,1) = zero
               vel(ulo(1):domlo(1)-1,j,k,2) = bc_vw(bcv)
               vel(ulo(1):domlo(1)-1,j,k,3) = bc_ww(bcv)
            end if
         end do
      end do
   endif

   if (nrgt .gt. 0) then
      do k=ulo(3),uhi(3)
         do j=ulo(2),uhi(2)
            bcv = bct_ihi(j,k,2)
            if (bct_ihi(j,k,1) == psw_) then
               vel(domhi(1)+1:uhi(1),j,k,1) = zero
               vel(domhi(1)+1:uhi(1),j,k,2) = bc_vw(bcv)
               vel(domhi(1)+1:uhi(1),j,k,3) = bc_ww(bcv)
            end if
         end do
      end do
   endif

   if (nbot .gt. 0) then
      do k=ulo(3),uhi(3)
         do i=ulo(1),uhi(1)
            bcv = bct_jlo(i,k,2)
            if (bct_jlo(i,k,1) == psw_)then
               vel(i,ulo(2):domlo(2)-1,k,1) = bc_uw(bcv)
               vel(i,ulo(2):domlo(2)-1,k,2) = zero
               vel(i,ulo(2):domlo(2)-1,k,3) = bc_ww(bcv)
            end if
         end do
      end do
   endif

   if (ntop .gt. 0) then
      do k=ulo(3),uhi(3)
         do i=ulo(1),uhi(1)
            bcv = bct_jhi(i,k,2)
            if (bct_jhi(i,k,1) == psw_)then
               vel(i,domhi(2)+1:uhi(2),k,1) = bc_uw(bcv)
               vel(i,domhi(2)+1:uhi(2),k,2) = zero
               vel(i,domhi(2)+1:uhi(2),k,3) = bc_ww(bcv)
            end if
         end do
      end do
   endif

   if (ndwn .gt. 0) then
      do j=ulo(2),uhi(2)
         do i=ulo(1),uhi(1)
            bcv = bct_klo(i,j,2)
            if (bct_klo(i,j,1) == psw_) then
               vel(i,j,ulo(3):domlo(3)-1,1) = bc_uw(bcv)
               vel(i,j,ulo(3):domlo(3)-1,2) = bc_vw(bcv)
               vel(i,j,ulo(3):domlo(3)-1,3) = zero
            end if
         end do
      end do
   endif

   if (nup .gt. 0) then
      do j=ulo(2),uhi(2)
         do i=ulo(1),uhi(1)
            bcv = bct_khi(i,j,2)
            if (bct_khi(i,j,1) == psw_) then
               vel(i,j,domhi(3)+1:uhi(3),1) = bc_uw(bcv)
               vel(i,j,domhi(3)+1:uhi(3),2) = bc_vw(bcv)
               vel(i,j,domhi(3)+1:uhi(3),3) = zero
            end if
         end do
      end do
   endif

end subroutine set_velocity_bcs

! *****************************************************************
