! 
!              
!  This subroutine sets the BCs for velocity components only.
!  
!  Author: Michele Rosso
! 
!  Date: December 20, 2017
!
! 
subroutine set_mac_velocity_bcs ( slo, shi, u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
     & bct_ilo, bct_ihi, bct_jlo, bct_jhi, bct_klo, bct_khi,               &
     & domlo, domhi, ng ) bind(C) 

   use amrex_fort_module,  only: ar => amrex_real
   use iso_c_binding ,     only: c_int
   use param,              only: zero, one, two, half
   use bc

   implicit none

   ! Array bounds
   integer(c_int), intent(in   ) :: slo(3), shi(3)
   integer(c_int), intent(in   ) :: ulo(3), uhi(3)
   integer(c_int), intent(in   ) :: vlo(3), vhi(3)
   integer(c_int), intent(in   ) :: wlo(3), whi(3)

   ! Grid bounds
   integer(c_int), intent(in   ) :: domlo(3), domhi(3)

   ! Number of ghost nodes
   integer(c_int), intent(in   ) :: ng
   
   ! BCs type
   integer(c_int), intent(in   ) :: &
        bct_ilo(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
        bct_ihi(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
        bct_jlo(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
        bct_jhi(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
        bct_klo(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
        bct_khi(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

   ! Arrays
   real(ar),      intent(inout) ::  &
        u_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)),      &
        v_g(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)),      &
        w_g(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

   ! Local variables
   integer  :: bcv, i, j, k
   integer  :: nlft, nrgt, nbot, ntop, nup, ndwn
   integer  :: num_cells_y
   real(ar) :: y

   nlft = max(0,domlo(1)-slo(1))
   nbot = max(0,domlo(2)-slo(2))
   ndwn = max(0,domlo(3)-slo(3))

   nrgt = max(0,shi(1)-domhi(1))
   ntop = max(0,shi(2)-domhi(2))
   nup  = max(0,shi(3)-domhi(3))

   num_cells_y = domhi(2) - domlo(2) + 1

   if (nlft .gt. 0) then
      do k = slo(3), shi(3)
         do j = slo(2), shi(2)
            y = (real(j,ar) + half) / num_cells_y

            bcv = bct_ilo(j,k,2)

            select case (bct_ilo(j,k,1))
               
            case ( pinf_, pout_) 
               
               u_g(ulo(1):domlo(1)-1,j,k) =      u_g(domlo(1),j,k)
               v_g(vlo(1):domlo(1)-1,j,k) =      v_g(domlo(1),j,k)
               w_g(wlo(1):domlo(1)-1,j,k) =      w_g(domlo(1),j,k)

            case ( minf_)


               u_g(ulo(1):domlo(1)  ,j,k) = 6 * bc_u_g(bcv) * y * (one - y)
               v_g(vlo(1):domlo(1)-1,j,k) = 0.0d0
               w_g(wlo(1):domlo(1)-1,j,k) = 0.0d0

            case ( nsw_) 

               u_g(ulo(1):domlo(1)  ,j,k) =  0.0d0
               v_g(vlo(1):domlo(1)-1,j,k) = -v_g(domlo(1),j,k)
               w_g(wlo(1):domlo(1)-1,j,k) = -w_g(domlo(1),j,k)

            case ( fsw_)
               
               u_g(ulo(1):domlo(1)  ,j,k) = 0.0d0
               v_g(vlo(1):domlo(1)-1,j,k) = v_g(domlo(1),j,k)
               w_g(wlo(1):domlo(1)-1,j,k) = w_g(domlo(1),j,k)

            end select
            
         end do
      end do
   endif

   if (nrgt .gt. 0) then
      
      do k = slo(3),shi(3)
         do j = slo(2),shi(2)
            
            bcv = bct_ihi(j,k,2)
            
            select case ( bct_ihi(j,k,1) )

            case ( pinf_, pout_ )
               
               u_g(domhi(1)+2:uhi(1),j,k) =      u_g(domhi(1)+1,j,k)
               v_g(domhi(1)+1:vhi(1),j,k) =      v_g(domhi(1)  ,j,k)
               w_g(domhi(1)+1:whi(1),j,k) =      w_g(domhi(1)  ,j,k)

            case ( minf_ )
               

               u_g(domhi(1)+1:uhi(1),j,k) = bc_u_g(bcv)
               v_g(domhi(1)+1:vhi(1),j,k) = 0.0d0
               w_g(domhi(1)+1:whi(1),j,k) = 0.0d0

            case ( nsw_ ) 

               u_g(domhi(1)+1:uhi(1),j,k) =  0.0d0
               v_g(domhi(1)+1:vhi(1),j,k) = -v_g(domhi(1),j,k)
               w_g(domhi(1)+1:whi(1),j,k) = -w_g(domhi(1),j,k)

            case ( fsw_ ) 

               u_g(domhi(1)+1:uhi(1),j,k) = 0.0d0
               v_g(domhi(1)+1:vhi(1),j,k) = v_g(domhi(1),j,k)
               w_g(domhi(1)+1:whi(1),j,k) = w_g(domhi(1),j,k)

            end select

         end do
      end do
   endif

   if (nbot .gt. 0) then
      
      do k = slo(3), shi(3)
         do i = slo(1), shi(1)
            
            bcv = bct_jlo(i,k,2)

            select case ( bct_jlo(i,k,1) )

            case ( pinf_, pout_) 

               u_g(i,ulo(2):domlo(2)-1,k) =      u_g(i,domlo(2),k)
               v_g(i,vlo(2):domlo(2)-1,k) =      v_g(i,domlo(2),k)
               w_g(i,wlo(2):domlo(2)-1,k) =      w_g(i,domlo(2),k)

            case ( minf_ )

               u_g(i,ulo(2):domlo(2)-1,k) = 0.0d0
               v_g(i,vlo(2):domlo(2)  ,k) = bc_v_g(bcv)
               w_g(i,wlo(2):domlo(2)-1,k) = 0.0d0

            case ( nsw_ )

               u_g(i,ulo(2):domlo(2)-1,k) = -u_g(i,domlo(2),k)
               v_g(i,vlo(2):domlo(2)  ,k) =  0.0d0
               w_g(i,wlo(2):domlo(2)-1,k) = -w_g(i,domlo(2),k)

            case ( fsw_)

               u_g(i,ulo(2):domlo(2)-1,k) = u_g(i,domlo(2),k)
               v_g(i,vlo(2):domlo(2)  ,k) = 0.0d0
               w_g(i,wlo(2):domlo(2)-1,k) = w_g(i,domlo(2),k)

            end select

         end do
      end do
   endif

   if (ntop .gt. 0) then

      do k = slo(3), shi(3)
         do i = slo(1), shi(1)
            
            bcv = bct_jhi(i,k,2)

            select case ( bct_jhi(i,k,1) )

            case ( pinf_, pout_ )
               
               u_g(i,domhi(2)+1:uhi(2),k) =      u_g(i,domhi(2)  ,k)
               v_g(i,domhi(2)+2:vhi(2),k) =      v_g(i,domhi(2)+1,k)
               w_g(i,domhi(2)+1:whi(2),k) =      w_g(i,domhi(2)  ,k)

            case ( minf_) 


               u_g(i,domhi(2)+1:uhi(2),k) = 0.0d0
               v_g(i,domhi(2)+1:vhi(2),k) = bc_v_g(bcv)
               w_g(i,domhi(2)+1:whi(2),k) = 0.0d0

            case ( nsw_) 

               u_g(i,domhi(2)+1:uhi(2),k) = -u_g(i,domhi(2),k)
               v_g(i,domhi(2)+1:vhi(2),k) =  0.0d0
               w_g(i,domhi(2)+1:whi(2),k) = -w_g(i,domhi(2),k)

            case ( fsw_)

               u_g(i,domhi(2)+1:uhi(2),k) = u_g(i,domhi(2),k)
               v_g(i,domhi(2)+1:vhi(2),k) = 0.0d0
               w_g(i,domhi(2)+1:whi(2),k) = w_g(i,domhi(2),k)

            end select
         end do
      end do
   endif

   if (ndwn .gt. 0) then

      do j = slo(2), shi(2)
         do i = slo(1), shi(1)
            
            bcv = bct_klo(i,j,2)

            select case (bct_klo(i,j,1))

            case ( pinf_, pout_ ) 

               u_g(i,j,ulo(3):domlo(3)-1) =      u_g(i,j,domlo(3))
               v_g(i,j,vlo(3):domlo(3)-1) =      v_g(i,j,domlo(3))
               w_g(i,j,wlo(3):domlo(3)-1) =      w_g(i,j,domlo(3))

            case ( minf_ )

               u_g(i,j,ulo(3):domlo(3)-1) = 0.0d0
               v_g(i,j,vlo(3):domlo(3)-1) = 0.0d0
               w_g(i,j,wlo(3):domlo(3)  ) = bc_w_g(bcv)

            case ( nsw_ )

               u_g(i,j,ulo(3):domlo(3)-1) = -u_g(i,j,domlo(3))
               v_g(i,j,vlo(3):domlo(3)-1) = -v_g(i,j,domlo(3))
               w_g(i,j,wlo(3):domlo(3)  ) =  0.0d0

            case ( fsw_ ) 

               u_g(i,j,ulo(3):domlo(3)-1) = u_g(i,j,domlo(3))
               v_g(i,j,vlo(3):domlo(3)-1) = v_g(i,j,domlo(3))
               w_g(i,j,wlo(3):domlo(3)  ) = 0.0d0
               
            end select
         end do
      end do
   endif

   if (nup .gt. 0) then

      do j = slo(2), shi(2)
         do i = slo(1), shi(1)
            
            bcv = bct_khi(i,j,2)

            select case ( bct_khi(i,j,1) )

            case ( pinf_, pout_ )
               
               u_g(i,j,domhi(3)+1:uhi(3)) =      u_g(i,j,domhi(3)  )
               v_g(i,j,domhi(3)+1:vhi(3)) =      v_g(i,j,domhi(3)  )
               w_g(i,j,domhi(3)+2:whi(3)) =      w_g(i,j,domhi(3)+1)
             
            case ( minf_ ) 

               u_g(i,j,domhi(3)+1:uhi(3)) = 0.0d0
               v_g(i,j,domhi(3)+1:vhi(3)) = 0.0d0
               w_g(i,j,domhi(3)+1:whi(3)) = bc_w_g(bcv)

            case ( nsw_ ) 

               u_g(i,j,domhi(3)+1:uhi(3)) = -u_g(i,j,domhi(3))
               v_g(i,j,domhi(3)+1:vhi(3)) = -v_g(i,j,domhi(3))
               w_g(i,j,domhi(3)+1:whi(3)) =  0.0d0

            case ( fsw_ ) 

               u_g(i,j,domhi(3)+1:uhi(3)) = u_g(i,j,domhi(3))
               v_g(i,j,domhi(3)+1:vhi(3)) = v_g(i,j,domhi(3))
               w_g(i,j,domhi(3)+1:whi(3)) = 0.0d0

            end select
         end do
      end do
   endif

   ! *********************************************************************************
   ! We have to do the PSW bc's last because otherwise non-zero moving wall values
   ! can get over-written
   ! *********************************************************************************

   if (nlft .gt. 0) then
      do k=slo(3),shi(3)
         do j=slo(2),shi(2)
            bcv = bct_ilo(j,k,2)
            if (bct_ilo(j,k,1) == PSW_) then
               u_g(ulo(1):domlo(1)  ,j,k) = 0.0d0
               v_g(vlo(1):domlo(1)-1,j,k) = 2.0*bc_vw_g(bcv) - v_g(domlo(1),j,k)
               w_g(wlo(1):domlo(1)-1,j,k) = 2.0*bc_ww_g(bcv) - w_g(domlo(1),j,k)
            end if
         end do
      end do
   endif

   if (nrgt .gt. 0) then
      do k=slo(3),shi(3)
         do j=slo(2),shi(2)
            bcv = bct_ihi(j,k,2)
            if (bct_ihi(j,k,1) == PSW_) then
               u_g(domhi(1)+1:uhi(1),j,k) = 0.0d0
               v_g(domhi(1)+1:vhi(1),j,k) = 2.0*bc_vw_g(bcv) - v_g(domhi(1),j,k)
               w_g(domhi(1)+1:whi(1),j,k) = 2.0*bc_ww_g(bcv) - w_g(domhi(1),j,k)
            end if
         end do
      end do
   endif

   if (nbot .gt. 0) then
      do k=slo(3),shi(3)
         do i=slo(1),shi(1)
            bcv = bct_jlo(i,k,2)
            if (bct_jlo(i,k,1) == PSW_)then
               u_g(i,ulo(2):domlo(2)-1,k) = 2.0*bc_uw_g(bcv) - u_g(i,domlo(2),k)
               v_g(i,vlo(2):domlo(2)  ,k) = 0.0d0
               w_g(i,wlo(2):domlo(2)-1,k) = 2.0*bc_ww_g(bcv) - w_g(i,domlo(2),k)
            end if
         end do
      end do
   endif

   if (ntop .gt. 0) then
      do k=slo(3),shi(3)
         do i=slo(1),shi(1)
            bcv = bct_jhi(i,k,2)
            if (bct_jhi(i,k,1) == PSW_)then
               u_g(i,domhi(2)+1:uhi(2),k) = 2.0*bc_uw_g(bcv) - u_g(i,domhi(2),k)
               v_g(i,domhi(2)+1:vhi(2),k) = 0.0d0
               w_g(i,domhi(2)+1:whi(2),k) = 2.0*bc_ww_g(bcv) - w_g(i,domhi(2),k)
            end if
         end do
      end do
   endif

   if (ndwn .gt. 0) then
      do j=slo(2),shi(2)
         do i=slo(1),shi(1)
            bcv = bct_klo(i,j,2)
            if (bct_klo(i,j,1) == PSW_) then
               u_g(i,j,ulo(3):domlo(3)-1) = 2.0*bc_uw_g(bcv) - u_g(i,j,domlo(3))
               v_g(i,j,vlo(3):domlo(3)-1) = 2.0*bc_vw_g(bcv) - v_g(i,j,domlo(3))
               w_g(i,j,wlo(3):domlo(3)  ) = 0.0d0
            end if
         end do
      end do
   endif

   if (nup .gt. 0) then
      do j=slo(2),shi(2)
         do i=slo(1),shi(1)
            bcv = bct_khi(i,j,2)
            if (bct_khi(i,j,1) == PSW_) then
               u_g(i,j,domhi(3)+1:uhi(3)) = 2.0*bc_uw_g(bcv) - u_g(i,j,domhi(3))
               v_g(i,j,domhi(3)+1:vhi(3)) = 2.0*bc_vw_g(bcv) - v_g(i,j,domhi(3))
               w_g(i,j,domhi(3)+1:whi(3)) = 0.0d0
            end if
         end do
      end do
   endif

end subroutine set_mac_velocity_bcs

