module eb_interpolation_mod

   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding,     only: c_int

   use constant,          only: zero, one

   implicit none
   private

   public interp_to_face_centroid

contains

   !
   ! Returns the flux at the face centroid ALREADY multiplied by the face area
   !
   function interp_to_face_centroid ( i, j, k, dir, var, vlo,  n,  &
                                     afrac, alo, cent, clo, nbr )  result(ivar)

      use amrex_ebcellflag_module, only: is_covered_cell
      use amrex_error_module,      only: amrex_abort

      ! Face indeces: these must be consistent with a staggered indexing
      ! and therefore consistent with the value of dir
      integer(c_int),  intent(in   ) :: i, j, k

      ! Direction of staggering (1=x, 2=y, 3=z): this specify how (i,j,k) must
      ! be interpreted, i.e. which staggered numbering the indexing refer to
      integer(c_int),  intent(in   ) :: dir

      ! The component to interpolate
      integer(c_int),  intent(in   ) :: n

      ! Array Bounds ( only start index )
      integer(c_int),  intent(in   ) :: vlo(3), alo(3), clo(3)

      ! Arrays
      real(ar),        intent(in   ) ::           &
           &   var(vlo(1):, vlo(2):, vlo(3):,1:), &
           & afrac(alo(1):, alo(2):, alo(3):),    &
           &  cent(clo(1):, clo(2):, clo(3):,1:)

      ! Neighbors information
      integer(c_int),  intent(in   ) :: nbr(-1:1,-1:1,-1:1)

      ! Output: the interpolated value
      real(ar)                       :: ivar

      ! Local variables
      real(ar)                       :: fracx, fracy, fracz

      select case ( dir )
      case(1) ! >>>>>>>>>>>>>>>>>>>>>>  X-face <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

         if ( afrac(i,j,k) == zero ) then
            ivar = zero
         else if ( afrac(i,j,k) == one ) then
            ivar = var(i,j,k,n)
         else
            if ( cent(i,j,k,1) < zero ) then
               fracy = - cent(i,j,k,1) * nbr(0,-1,0)
               if ( cent(i,j,k,2) <= zero ) then
                  fracz = - cent(i,j,k,2) * nbr(0,0,-1)
                  ivar = (one-fracz) * (     fracy * var(i,j-1,k  ,n)  + &
                       &               (one-fracy) * var(i,j  ,k  ,n)) + &
                       &       fracz * (     fracy * var(i,j-1,k-1,n)  + &
                       &               (one-fracy) * var(i,j  ,k-1,n))
               else
                  fracz =  cent(i,j,k,2) * nbr(0,0,1)
                  ivar = (one-fracz) * (     fracy * var(i,j-1,k  ,n)  + &
                       &               (one-fracy) * var(i,j  ,k  ,n)) + &
                       &       fracz * (     fracy * var(i,j-1,k+1,n)  + &
                       &               (one-fracy) * var(i,j  ,k+1,n))
               endif
            else
               fracy = cent(i,j,k,1) * nbr(0,1,0)
               if ( cent(i,j,k,2) <= zero ) then
                  fracz = - cent(i,j,k,2) * nbr(0,0,-1)
                  ivar = (one-fracz) * (     fracy * var(i,j+1,k  ,n)  + &
                       &               (one-fracy) * var(i,j  ,k  ,n)) + &
                       &       fracz * (     fracy * var(i,j+1,k-1,n)  + &
                       &               (one-fracy) * var(i,j  ,k-1,n))
               else
                  fracz =  cent(i,j,k,2) * nbr(0,0,1)
                  ivar= (one-fracz) * (     fracy * var(i,j+1,k  ,n)  + &
                       &              (one-fracy) * var(i,j  ,k  ,n)) + &
                       &      fracz * (     fracy * var(i,j+1,k+1,n)  + &
                       &              (one-fracy) * var(i,j  ,k+1,n))
               endif
            end if
         end if

      case(2)  ! >>>>>>>>>>>>>>>>>>>>>>  Y-face <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

         if ( afrac(i,j,k) == zero ) then
            ivar = zero
         else if ( afrac(i,j,k) == one ) then
            ivar = var(i,j,k,n)
         else
            if ( cent(i,j,k,1) < zero ) then
               fracx = - cent(i,j,k,1) * nbr(-1,0,0)
               if ( cent(i,j,k,2) <= zero ) then
                  fracz = - cent(i,j,k,2) * nbr(0,0,-1)
                  ivar = (one-fracz) * (     fracx * var(i-1,j,k  ,n)  + &
                       &               (one-fracx) * var(i  ,j,k  ,n)) + &
                       &       fracz * (     fracx * var(i-1,j,k-1,n)  + &
                       &               (one-fracx) * var(i  ,j,k-1,n))
               else
                  fracz =  cent(i,j,k,2) * nbr(0,0,1)
                  ivar = (one-fracz) * (     fracx * var(i-1,j,k  ,n)  + &
                       &               (one-fracx) * var(i  ,j,k  ,n)) + &
                       &       fracz * (     fracx * var(i-1,j,k+1,n)  + &
                       &               (one-fracx) * var(i  ,j,k+1,n))
               endif
            else
               fracx = cent(i,j,k,1) * nbr(1,0,0)
               if ( cent(i,j,k,2) <= zero ) then
                  fracz = - cent(i,j,k,2) * nbr(0,0,-1)
                  ivar = (one-fracz) * (     fracx * var(i+1,j,k  ,n)  + &
                       &               (one-fracx) * var(i  ,j,k  ,n)) + &
                       &       fracz * (     fracx * var(i+1,j,k-1,n)  + &
                       &               (one-fracx) * var(i  ,j,k-1,n))
               else
                  fracz =  cent(i,j,k,2) * nbr(0,0,1)
                  ivar = (one-fracz) * (     fracx * var(i+1,j,k  ,n)  + &
                       &               (one-fracx) * var(i  ,j,k  ,n)) + &
                       &       fracz * (     fracx * var(i+1,j,k+1,n)  + &
                       &               (one-fracx) * var(i  ,j,k+1,n))
               endif
            end if
         end if

      case(3) ! >>>>>>>>>>>>>>>>>>>>>>  Z-face <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

         if ( afrac(i,j,k) == zero ) then
            ivar = zero
         else if ( afrac(i,j,k) == one ) then
            ivar = var(i,j,k,n)
         else
            if ( cent(i,j,k,1) < zero ) then
               fracx = - cent(i,j,k,1) * nbr(-1,0,0)
               if ( cent(i,j,k,2) <= zero ) then
                  fracy = - cent(i,j,k,2) * nbr(0,-1,0)
                  ivar = (one-fracy) * (     fracx * var(i-1,j  ,k,n)  + &
                       &               (one-fracx) * var(i  ,j  ,k,n)) + &
                       &       fracy * (     fracx * var(i-1,j-1,k,n)  + &
                       &               (one-fracx) * var(i  ,j-1,k,n))
               else
                  fracy =  cent(i,j,k,2) * nbr(0,1,0)
                  ivar = (one-fracy) * (     fracx * var(i-1,j  ,k,n)  + &
                       &               (one-fracx) * var(i  ,j  ,k,n)) + &
                       &       fracy * (     fracx * var(i-1,j+1,k,n)  + &
                       &               (one-fracx) * var(i  ,j+1,k,n))
               endif
            else
               fracx = cent(i,j,k,1) * nbr(1,0,0)
               if ( cent(i,j,k,2) <= zero ) then
                  fracy = - cent(i,j,k,2) * nbr(0,-1,0)
                  ivar = (one-fracy) * (     fracx * var(i+1,j  ,k,n)  + &
                       &               (one-fracx) * var(i  ,j  ,k,n)) + &
                       &       fracy * (     fracx * var(i+1,j-1,k,n)  + &
                       &               (one-fracx) * var(i  ,j-1,k,n))
               else
                  fracy =  cent(i,j,k,2) * nbr(0,1,0)
                  ivar = (one-fracy) * (     fracx * var(i+1,j  ,k,n)  + &
                       &               (one-fracx) * var(i  ,j  ,k,n)) + &
                       &       fracy * (     fracx * var(i+1,j+1,k,n)  + &
                       &               (one-fracx) * var(i  ,j+1,k,n))
               endif
            end if
         end if

      case default

         call amrex_abort( "interp_to_face_centroid(): value of 'dir'"&
                          //" is invalid. Must be either 1,2, or 3")

      end select

      !
      ! Return the flux multiplied by the face area
      !
      ivar = ivar * afrac(i,j,k)

   end function interp_to_face_centroid

end module eb_interpolation_mod
