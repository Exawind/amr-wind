! 
!              
!  This subroutine sets the BCs for all scalar variables involved in
!  the projection, EXCEPT the pressure and velocity components. 
!  
!  Author: Michele Rosso
! 
!  Date: December 20, 2017
!
! 
subroutine set_scalar_bcs ( ro, mu, lambda, slo, shi, &
     & bct_ilo, bct_ihi, bct_jlo, bct_jhi, bct_klo, bct_khi,           &
     & domlo, domhi, ng ) bind(C) 

   use amrex_fort_module,  only: ar => amrex_real
   use iso_c_binding ,     only: c_int
   use bc
   use constant,           only: ro_0, mu_0
   use param    ,          only: is_undefined
   
   implicit none

   ! Array bounds
   integer(c_int), intent(in   ) :: slo(3), shi(3)

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
        ro(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),     &
        mu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),     &
        lambda(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

   ! Local variables
   integer  :: bcv, i, j, k
   integer  :: nlft, nrgt, nbot, ntop, nup, ndwn
   real(ar) :: bc_ro, bc_mu, bc_lambda

   nlft = max(0,domlo(1)-slo(1))
   nbot = max(0,domlo(2)-slo(2))
   ndwn = max(0,domlo(3)-slo(3))

   nrgt = max(0,shi(1)-domhi(1))
   ntop = max(0,shi(2)-domhi(2))
   nup  = max(0,shi(3)-domhi(3))

   bc_ro = ro_0
   
   if (nlft .gt. 0) then
      do k = slo(3), shi(3)
         do j = slo(2), shi(2)

            bcv = bct_ilo(j,k,2)

            select case (bct_ilo(j,k,1))
               
            case ( pinf_, pout_, nsw_, fsw_) 

               ro(slo(1):domlo(1)-1,j,k) =     ro(domlo(1),j,k)
               mu(slo(1):domlo(1)-1,j,k) =     mu(domlo(1),j,k)
               lambda(slo(1):domlo(1)-1,j,k) = lambda(domlo(1),j,k)

            case ( minf_)

               bc_mu     = mu_0
               bc_lambda = -(2.0d0/3.0d0) * mu_0

               ro(slo(1):domlo(1)-1,j,k)     = bc_ro
               mu(slo(1):domlo(1)-1,j,k)     = bc_mu
               lambda(slo(1):domlo(1)-1,j,k) = bc_lambda

            end select
            
         end do
      end do
   endif

   if (nrgt .gt. 0) then
      
      do k = slo(3),shi(3)
         do j = slo(2),shi(2)
            
            bcv = bct_ihi(j,k,2)
            
            select case ( bct_ihi(j,k,1) )

            case ( pinf_, pout_, nsw_, fsw_) 

                   ro(domhi(1)+1:shi(1),j,k) =     ro(domhi(1)  ,j,k)
                   mu(domhi(1)+1:shi(1),j,k) =     mu(domhi(1)  ,j,k)
               lambda(domhi(1)+1:shi(1),j,k) = lambda(domhi(1)  ,j,k)

            case ( minf_ )

               bc_mu     = mu_0
               bc_lambda = -(2.0d0/3.0d0) * mu_0
               
                   ro(domhi(1)+1:shi(1),j,k) = bc_ro
                   mu(domhi(1)+1:shi(1),j,k) = bc_mu
               lambda(domhi(1)+1:shi(1),j,k) = bc_lambda

            end select

         end do
      end do
   endif

   if (nbot .gt. 0) then
      
      do k = slo(3), shi(3)
         do i = slo(1), shi(1)
            
            bcv = bct_jlo(i,k,2)

            select case ( bct_jlo(i,k,1) )

            case ( pinf_, pout_, nsw_, fsw_) 

                   ro(i,slo(2):domlo(2)-1,k) =     ro(i,domlo(2),k)
                   mu(i,slo(2):domlo(2)-1,k) =     mu(i,domlo(2),k)
               lambda(i,slo(2):domlo(2)-1,k) = lambda(i,domlo(2),k)

            case ( minf_ )

               bc_mu     = mu_0
               bc_lambda = -(2.0d0/3.0d0) * mu_0

                   ro(i,slo(2):domlo(2)-1,k) = bc_ro
                   mu(i,slo(2):domlo(2)-1,k) = bc_mu
               lambda(i,slo(2):domlo(2)-1,k) = bc_lambda

            end select

         end do
      end do
   endif

   if (ntop .gt. 0) then

      do k = slo(3), shi(3)
         do i = slo(1), shi(1)
            
            bcv = bct_jhi(i,k,2)

            select case ( bct_jhi(i,k,1) )

            case ( pinf_, pout_, nsw_, fsw_) 

                   ro(i,domhi(2)+1:shi(2),k) =     ro(i,domhi(2)  ,k)
                   mu(i,domhi(2)+1:shi(2),k) =     mu(i,domhi(2)  ,k)
               lambda(i,domhi(2)+1:shi(2),k) = lambda(i,domhi(2)  ,k)

            case ( minf_) 

               bc_mu     = mu_0
               bc_lambda = -(2.0d0/3.0d0) * mu_0

                   ro(i,domhi(2)+1:shi(2),k) = bc_ro
                   mu(i,domhi(2)+1:shi(2),k) = bc_mu
               lambda(i,domhi(2)+1:shi(2),k) = bc_lambda

            end select
         end do
      end do
   endif

   if (ndwn .gt. 0) then

      do j = slo(2), shi(2)
         do i = slo(1), shi(1)
            
            bcv = bct_klo(i,j,2)

            select case (bct_klo(i,j,1))

            case ( pinf_, pout_, nsw_, fsw_) 

                   ro(i,j,slo(3):domlo(3)-1) =     ro(i,j,domlo(3))
                   mu(i,j,slo(3):domlo(3)-1) =     mu(i,j,domlo(3))
               lambda(i,j,slo(3):domlo(3)-1) = lambda(i,j,domlo(3))

            case ( minf_ )

               bc_mu     = mu_0
               bc_lambda = -(2.0d0/3.0d0) * mu_0
               
                   ro(i,j,slo(3):domlo(3)-1) = bc_ro
                   mu(i,j,slo(3):domlo(3)-1) = bc_mu
               lambda(i,j,slo(3):domlo(3)-1) = bc_lambda

            end select
         end do
      end do
   endif

   if (nup .gt. 0) then

      do j = slo(2), shi(2)
         do i = slo(1), shi(1)
            
            bcv = bct_khi(i,j,2)

            select case ( bct_khi(i,j,1) )

            case ( pinf_, pout_, nsw_, fsw_) 
             
                   ro(i,j,domhi(3)+1:shi(3)) =     ro(i,j,domhi(3)  )
                   mu(i,j,domhi(3)+1:shi(3)) =     mu(i,j,domhi(3)  )
               lambda(i,j,domhi(3)+1:shi(3)) = lambda(i,j,domhi(3)  )

            case ( minf_ ) 

               bc_mu     = mu_0
               bc_lambda = -(2.0d0/3.0d0) * mu_0
               
                   ro(i,j,domhi(3)+1:shi(3)) = bc_ro
                   mu(i,j,domhi(3)+1:shi(3)) = bc_mu
               lambda(i,j,domhi(3)+1:shi(3)) = bc_lambda

            end select
         end do
      end do
   endif

end subroutine set_scalar_bcs
