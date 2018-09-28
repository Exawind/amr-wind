! 
!              
!  This subroutine sets the BCs for the pressure field.
! 
!  Most of the times, this routines will simply perform
!  0th-order extrapolation of the pressure values from
!  iterior cells to exterior cells. This corresponds to
!  Neumann's BCs on the pressure. 
!
!  When the pressure is periodic in a certain direction,
!  nothing will be done here along that same direction.
!  Periodicity is handled on the C++ side.
!
!  For all the boundary conditions demanding a Dirichlet's
!  boundary for pressure, this routine will perform
!  1st-order extrapolation of the pressure value from the
!  interior cells to the exterior cells.
!  In particular, when the pressure is imposed at an inlet
!  /outlet, the pressure at the first external cell is
!  extrapolated by using the given pressure boundary value
!  and the pressure value at the first internal cell.
!  Instead, when a boundary pressure value must be computed
!  must be computed a part of a mass-inflow condition,
!  the extrapolation is performed by using the pressure
!  at the first two interior cells.
!  
!  Author: Michele Rosso
! 
!  Date: December 20, 2017
!
! 
subroutine extrap_pressure_to_ghost_cells (  phi, slo, shi, bct_ilo, bct_ihi, &
     bct_jlo, bct_jhi, bct_klo, bct_khi, domlo, domhi, ng ) bind(C) 

   use amrex_fort_module,  only: ar => amrex_real
   use iso_c_binding ,     only: c_int
   use bc
   
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
        phi(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))


   ! Local variables
   integer             :: bcv, i, j, k
   integer             :: nlft, nrgt, nbot, ntop, nup, ndwn
   real(ar), parameter :: two = 2.0_ar 

   nlft = max(0,domlo(1)-slo(1))
   nbot = max(0,domlo(2)-slo(2))
   ndwn = max(0,domlo(3)-slo(3))

   nrgt = max(0,shi(1)-domhi(1))
   ntop = max(0,shi(2)-domhi(2))
   nup  = max(0,shi(3)-domhi(3))


   !  X low
   if ( domlo(1) > slo(1) .and. .not. cyclic_x ) then
      
      do k = slo(3), shi(3)
         do j = slo(2), shi(2)

            bcv = bct_ilo(j,k,2)

            select case ( bct_ilo(j,k,1) )
               
            case ( pinf_, pout_) 
     
               ! Dirichlet value is set to face: extrapolate value at
               ! first ghost cell
               phi(domlo(1)-1,j,k) =  2.d0*phi(domlo(1)-1,j,k) - phi(domlo(1),j,k)

            case default

               phi(slo(1):domlo(1)-1,j,k) = phi(domlo(1),j,k)

            end select
            
         end do
      end do
   endif
   
   ! X high
   if ( domhi(1) < shi(1) .and. .not. cyclic_x ) then
      
      do k = slo(3),shi(3)
         do j = slo(2),shi(2)
            
            bcv = bct_ihi(j,k,2)
            
            select case ( bct_ihi(j,k,1) )

            case ( pinf_, pout_ )

               ! Dirichlet value is set to face: extrapolate value at
               ! first ghost cell
               phi(domhi(1)+1,j,k) = 2.d0*phi(domhi(1)+1,j,k) - phi(domhi(1),j,k)

            case default
               
               phi(domhi(1)+1:shi(1),j,k) =  phi(domhi(1),j,k) 

            end select

         end do
      end do
   endif

   ! Y low
   if ( domlo(2) > slo(2) .and. .not. cyclic_y ) then
      
      do k = slo(3), shi(3)
         do i = slo(1), shi(1)
            
            bcv = bct_jlo(i,k,2)

            select case ( bct_jlo(i,k,1) )

            case ( pinf_, pout_) 

               ! Dirichlet value is set to face: extrapolate value at
               ! first ghost cell
               phi(i,domlo(2)-1,k) =  2.d0*phi(i,domlo(2)-1,k) - phi(i,domlo(2),k)

            case default 

               phi(i,slo(2):domlo(2)-1,k) = phi(i,domlo(2),k)

            end select

         end do
      end do
   endif

   ! Y high
   if ( domhi(2) < shi(2) .and. .not. cyclic_y ) then

      do k = slo(3), shi(3)
         do i = slo(1), shi(1)
            
            bcv = bct_jhi(i,k,2)

            select case ( bct_jhi(i,k,1) )

            case ( pinf_, pout_ )

               ! Dirichlet value is set to face: extrapolate value at
               ! first ghost cell
               phi(i,domhi(2)+1,k) =  2.d0*phi(i,domhi(2)+1,k) - phi(i,domhi(2),k)

            case default 

               phi(i,domhi(2)+1:shi(2),k) = phi(i,domhi(2),k) 

            end select

         end do
      end do
   endif


   ! Z low
   if ( domlo(3) > slo(3) .and. .not. cyclic_z ) then

      do j = slo(2), shi(2)
         do i = slo(1), shi(1)
            
            bcv = bct_klo(i,j,2)

            select case (bct_klo(i,j,1))

            case ( pinf_, pout_ ) 

               ! Dirichlet value is set to face: extrapolate value at
               ! first ghost cell
               phi(i,j,domlo(3)-1) =  2.d0*phi(i,j,domlo(3)-1) - phi(i,j,domlo(3))

            case default

               phi(i,j,slo(3):domlo(3)-1) = phi(i,j,domlo(3))
               
            end select
            
         end do
      end do
   endif


   ! Z high
   if ( domhi(3) < shi(3) .and. .not. cyclic_z ) then

      do j = slo(2), shi(2)
         do i = slo(1), shi(1)
            
            bcv = bct_khi(i,j,2)

            select case ( bct_khi(i,j,1) )

            case ( pinf_, pout_ )
               
               ! Dirichlet value is set to face: extrapolate value at
               ! first ghost cell
               phi(i,j,domhi(3)+1) =  2.d0*phi(i,j,domhi(3)+1) - phi(i,j,domhi(3))

            case default 

               phi(i,j,domhi(3)+1:shi(3)) = phi(i,j,domhi(3)) 

            end select
            
         end do
      end do
   endif

end subroutine extrap_pressure_to_ghost_cells
