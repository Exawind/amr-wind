module eb_interpolation_mod

   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one
   use bc,                only: minf_, nsw_, fsw_, psw_, pinf_, pout_

   implicit none
   private

   public interpolate_to_face_centroid

contains


   ! Interpolate x-face variable from face center to face centroid  
   subroutine interpolate_to_face_centroid ( lo, hi, ivar, var, vlo, vhi, ncomp, &
        areafrac, alo, ahi, cent, clo, chi, flags, flo, fhi, face_type  ) 

      use amrex_ebcellflag_module, only: is_covered_cell, get_neighbor_cells

      ! Tile bounds ( face centered )
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array Bounds
      integer(c_int),  intent(in   ) :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: alo(3), ahi(3)
      integer(c_int),  intent(in   ) :: clo(3), chi(3)
      integer(c_int),  intent(in   ) :: flo(3), fhi(3)
      integer(c_int),  intent(in   ) :: ncomp

      ! Type of face (1=x, 2=y, 3=z)
      integer(c_int),  intent(in   ) :: face_type


      ! Arrays
      real(ar),        intent(inout) ::                            &
           & ivar(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),ncomp)  ! Interpolated Variable

      real(ar),        intent(in   ) ::                            &
           & var(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),ncomp), &
           & areafrac(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3)),  &
           & cent(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3),2)

      integer(c_int),  intent(in   ) :: &
           & flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))


      ! Local variables
      integer(c_int)                 :: i, j, k, n, nbr(-1:1,-1:1,-1:1)
      real(ar)                       :: fracx, fracy, fracz

      select case ( face_type )
      case(1) ! >>>>>>>>>>>>>>>>>>>>>>  X-face <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !
         do n = 1, ncomp
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if ( ( areafrac(i,j,k) > zero ) .and. ( areafrac(i,j,k) < one ) ) then

                        call get_neighbor_cells( flags(i,j,k), nbr )

                        if ( cent(i,j,k,1) < zero ) then
                           fracy = - cent(i,j,k,1) * nbr(0,-1,0)                      
                           if ( cent(i,j,k,2) <= zero ) then
                              fracz = - cent(i,j,k,2) * nbr(0,0,-1)
                              ivar(i,j,k,n) = (one-fracz) * (     fracy * var(i,j-1,k  ,n)  + &
                                   &                        (one-fracy) * var(i,j  ,k  ,n)) + &
                                   &                fracz * (     fracy * var(i,j-1,k-1,n)  + &
                                   &                        (one-fracy) * var(i,j  ,k-1,n))
                           else
                              fracz =  cent(i,j,k,2) * nbr(0,0,1)
                              ivar(i,j,k,n) = (one-fracz) * (     fracy * var(i,j-1,k  ,n)  + &
                                   &                        (one-fracy) * var(i,j  ,k  ,n)) + &
                                   &                fracz * (     fracy * var(i,j-1,k+1,n)  + &
                                   &                        (one-fracy) * var(i,j  ,k+1,n))
                           endif
                        else
                           fracy = cent(i,j,k,1) * nbr(0,1,0)
                           if ( cent(i,j,k,2) <= zero ) then
                              fracz = - cent(i,j,k,2) * nbr(0,0,-1)
                              ivar(i,j,k,n) = (one-fracz) * (     fracy * var(i,j+1,k  ,n)  + &
                                   &                        (one-fracy) * var(i,j  ,k  ,n)) + &
                                   &                fracz * (     fracy * var(i,j+1,k-1,n)  + &
                                   &                        (one-fracy) * var(i,j  ,k-1,n))
                           else
                              fracz =  cent(i,j,k,2) * nbr(0,0,1)
                              ivar(i,j,k,n) = (one-fracz) * (     fracy * var(i,j+1,k  ,n)  + &
                                   &                        (one-fracy) * var(i,j  ,k  ,n)) + &
                                   &                fracz * (     fracy * var(i,j+1,k+1,n)  + &
                                   &                        (one-fracy) * var(i,j  ,k+1,n)) 
                           endif
                        end if
                     else                  
                        ivar(i,j,k,n) = var(i,j,k,n)
                     end if
                  end do
               end do
            end do
         end do

      case(2)  ! >>>>>>>>>>>>>>>>>>>>>>  Y-face <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

         do n = 1, ncomp
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if ( ( areafrac(i,j,k) > zero ) .and. ( areafrac(i,j,k) < one ) ) then

                        call get_neighbor_cells( flags(i,j,k), nbr )

                        if ( cent(i,j,k,1) < zero ) then
                           fracx = - cent(i,j,k,1) * nbr(-1,0,0)                      
                           if ( cent(i,j,k,2) <= zero ) then
                              fracz = - cent(i,j,k,2) * nbr(0,0,-1)
                              ivar(i,j,k,n) = (one-fracz) * (     fracx * var(i-1,j,k  ,n)  + &
                                   &                        (one-fracx) * var(i-1,j,k  ,n)) + &
                                   &                fracz * (     fracx * var(i  ,j,k-1,n)  + &
                                   &                        (one-fracx) * var(i  ,j,k-1,n))
                           else
                              fracz =  cent(i,j,k,2) * nbr(0,0,1)
                              ivar(i,j,k,n) = (one-fracz) * (     fracx * var(i-1,j,k  ,n)  + &
                                   &                        (one-fracx) * var(i  ,j,k  ,n)) + &
                                   &                fracz * (     fracx * var(i-1,j,k+1,n)  + &
                                   &                        (one-fracx) * var(i  ,j,k+1,n))
                           endif
                        else
                           fracx = cent(i,j,k,1) * nbr(1,0,0)
                           if ( cent(i,j,k,2) <= zero ) then
                              fracz = - cent(i,j,k,2) * nbr(0,0,-1)
                              ivar(i,j,k,n) = (one-fracz) * (     fracx * var(i+1,j,k  ,n)  + &
                                   &                        (one-fracx) * var(i  ,j,k  ,n)) + &
                                   &                fracz * (     fracx * var(i+1,j,k-1,n)  + &
                                   &                        (one-fracx) * var(i  ,j,k-1,n))
                           else
                              fracz =  cent(i,j,k,2) * nbr(0,0,1)
                              ivar(i,j,k,n) = (one-fracz) * (     fracx * var(i+1,j,k  ,n)  + &
                                   &                        (one-fracx) * var(i  ,j,k  ,n)) + &
                                   &                fracz * (     fracx * var(i+1,j,k+1,n)  + &
                                   &                        (one-fracx) * var(i  ,j,k+1,n)) 
                           endif
                        end if
                     else                  
                        ivar(i,j,k,n) = var(i,j,k,n)
                     end if
                  end do
               end do
            end do
         end do

      case(3) ! >>>>>>>>>>>>>>>>>>>>>>  Z-face <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !
         
         do n = 1, ncomp
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if ( ( areafrac(i,j,k) > zero ) .and. ( areafrac(i,j,k) < one ) ) then

                        call get_neighbor_cells( flags(i,j,k), nbr )

                        if ( cent(i,j,k,1) < zero ) then
                           fracx = - cent(i,j,k,1) * nbr(-1,0,0)                      
                           if ( cent(i,j,k,2) <= zero ) then
                              fracy = - cent(i,j,k,2) * nbr(0,-1,0)
                              ivar(i,j,k,n) = (one-fracy) * (     fracx * var(i-1,j  ,k,n)  + &
                                   &                        (one-fracx) * var(i  ,j  ,k,n)) + &
                                   &                fracy * (     fracx * var(i-1,j-1,k,n)  + &
                                   &                        (one-fracx) * var(i  ,j-1,k,n))
                           else
                              fracy =  cent(i,j,k,2) * nbr(0,1,0)
                              ivar(i,j,k,n) = (one-fracy) * (     fracx * var(i-1,j  ,k,n)  + &
                                   &                        (one-fracx) * var(i  ,j  ,k,n)) + &
                                   &                fracy * (     fracx * var(i-1,j+1,k,n)  + &
                                   &                        (one-fracx) * var(i  ,j+1,k,n))
                           endif
                        else
                           fracx = cent(i,j,k,1) * nbr(1,0,0)
                           if ( cent(i,j,k,2) <= zero ) then
                              fracy = - cent(i,j,k,2) * nbr(0,-1,0)
                              ivar(i,j,k,n) = (one-fracy) * (     fracx * var(i+1,j  ,k,n)  + &
                                   &                        (one-fracx) * var(i  ,j  ,k,n)) + &
                                   &                fracy * (     fracx * var(i+1,j-1,k,n)  + &
                                   &                        (one-fracx) * var(i  ,j-1,k,n))
                           else
                              fracy =  cent(i,j,k,2) * nbr(0,1,0)
                              ivar(i,j,k,n) = (one-fracy) * (     fracx * var(i+1,j  ,k,n)  + &
                                   &                        (one-fracx) * var(i  ,j  ,k,n)) + &
                                   &                fracy * (     fracx * var(i+1,j+1,k,n)  + &
                                   &                        (one-fracx) * var(i  ,j+1,k,n)) 
                           endif
                        end if
                     else                  
                        ivar(i,j,k,n) = var(i,j,k,n)
                     end if
                  end do
               end do
            end do
         end do
         
      case default
         
         write(*,*) "interpolate_to_face_centroid(): face_type = ", face_type, " but valid values are 1,2,3"
         stop
         
      end select

   end subroutine interpolate_to_face_centroid


end module eb_interpolation_mod
