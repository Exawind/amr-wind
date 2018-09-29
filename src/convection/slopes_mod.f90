!
!  This module contains the subroutines to compute the slopes in the
!  three directions for each velocity component.
!  The x,y, and z slopes for each velocity component are calculated at the
!  velocity component location via the second order Monotonized Central (MC)
!  limiter (van Leer, 1977). The scheme is described below for the u-velocity.
!
!
!                      |--x--|--x--|--x--|
!                        i-1    i    i+1
!
!  In the sketch above, the x represents the u-velocities while the vertical
!  bars | | enclose a u-cell (NOT a scalar cell!).
!  The MC limiter computes the slope at cell "i" by combining the left, central
!  and right u-variation "du":
!
!       du_l = u(i) - u(i-1)               = left variation
!       du_c = 0.5 * ( u(i+1) - u(i-1) )   = central (umlimited) variation
!       du_r = u(i+1) - u(i)               = right variation
!
!  Finally, the u-variation at cell "i" is given by :
!
!       du(i) = sign(du_c) * min(2|du_l|, |du_c|, 2|du_r|)) if du_l*du_r > 0
!       du(i) = 0                                           otherwise
!
!  The above procedure is applied direction by direction.
!
!  BOUNDARY CONDITIONS
!  When periodic or Neumann's BCs are imposed, the scheme can be applied
!  without any change since the ghost cells at the boundary are filled
!  by either periodicity or by extrapolation.
!  For Dirichlet's BCs in the transversal direction, the scheme can again
!  be applied as is since the velocity is known at the first ghost cell
!  out of the domain.
!  However, for Dirichlet's BCs in the longitudinal direction, the velocity
!  is not known outside the domain since the BC is applied directly at the first
!  valid node which lies on the boundary itself. Therefore, the scheme must be
!  arranged as follows to use ONLY values from inside the domain.
!  For a left boundary (i=0), the u-variations are:
!
!       du_l = 0                             Don't use values on the left
!       du_c = -1.5*u(0) + 2*u(1) -0.5*u(2)  2nd order right-biased
!       du_r = u(1) - u(0)                   Right variation
!
module slopes_mod

   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one
   use bc,                only: minf_, nsw_, fsw_, psw_

   implicit none
   private

   public compute_slopes

contains

   !
   ! Compute slopes
   !
   subroutine compute_slopes ( lo, hi, vel, vlo, vhi, &
                              xslopes, yslopes, zslopes, slo, shi, &
                              domlo, domhi, &
                              bc_ilo_type, bc_ihi_type, &
                              bc_jlo_type, bc_jhi_type, &
                              bc_klo_type, bc_khi_type, ng ) bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) :: lo(3), hi(3)
      integer(c_int), intent(in   ) :: ng

      ! Array bounds
      integer(c_int), intent(in   ) :: vlo(3), vhi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)

      ! Grid bounds
      integer(c_int), intent(in   ) :: domlo(3), domhi(3)

      ! BC types
      integer(c_int), intent(in   ) ::  &
           & bc_ilo_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_ihi_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_jlo_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_jhi_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_klo_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
           & bc_khi_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

      ! Arrays
      real(ar),       intent(in   ) ::                      &
           & vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3)

      real(ar),       intent(  out) ::                             &
           & xslopes(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3), &
           & yslopes(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3), &
           & zslopes(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      ! Local variables
      integer                       :: i, j, k, n
      real(ar)                      :: du_l, du_c, du_r
      real(ar),    parameter        :: two = 2.0_ar, three2nds = 1.5_ar

      do n = 1, 3
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! X direction
               du_l = two*(vel(i  ,j,k,n) - vel(i-1,j,k,n))
               du_r = two*(vel(i+1,j,k,n) - vel(i  ,j,k,n))
               du_c = half * ( vel(i+1,j,k,n) - vel(i-1,j,k,n) )
               xslopes(i,j,k,n) = mc_limiter ( du_l, du_c, du_r )

               ! Y direction
               du_l = two*(vel(i,j  ,k,n) - vel(i,j-1,k,n))
               du_r = two*(vel(i,j+1,k,n) - vel(i,j  ,k,n))
               du_c = half * ( vel(i,j+1,k,n) - vel(i,j-1,k,n) )
               yslopes(i,j,k,n) = mc_limiter ( du_l, du_c, du_r )

               ! z direction
               du_l = two*(vel(i,j,k  ,n) - vel(i,j,k-1,n))
               du_r = two*(vel(i,j,k+1,n) - vel(i,j,k  ,n))
               du_c = half * ( vel(i,j,k+1,n) - vel(i,j,k-1,n) )
               zslopes(i,j,k,n) = mc_limiter ( du_l, du_c, du_r )

            end do
         end do
      end do
      end do

      !
      ! Compute slopes at boundary where physical BCs are imposed
      !
      if ( lo(1) == domlo(1) ) then

         i = lo(1)

         do n = 1, 3
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)

               if ( ( bc_ilo_type(j,k,1) == MINF_ ) .or. &
                   ( bc_ilo_type(j,k,1) == NSW_ )  .or. &
                   ( bc_ilo_type(j,k,1) == FSW_ )  .or. &
                   ( bc_ilo_type(j,k,1) == PSW_ )  ) then

                  du_l = two*(vel(i  ,j,k,n) - vel(i-1,j,k,n))
                  du_r = two*(vel(i+1,j,k,n) - vel(i  ,j,k,n))
                  du_c = (vel(i+1,j,k,n)+3.d0*vel(i,j,k,n)-4.d0*vel(i-1,j,k,n))/3.d0
                  xslopes(i,j,k,n) = mc_limiter ( du_l, du_c, du_r )

               end if
            end do
         end do
         end do
      end if

      if ( hi(1) == domhi(1) ) then

         i = hi(1)

         do n = 1, 3
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)

               if ( ( bc_ihi_type(j,k,1) == MINF_ ) .or. &
                   ( bc_ihi_type(j,k,1) == NSW_ )  .or. &
                   ( bc_ihi_type(j,k,1) == FSW_ )  .or. &
                   ( bc_ihi_type(j,k,1) == PSW_ )  ) then

                  du_l = two*(vel(i  ,j,k,n) - vel(i-1,j,k,n))
                  du_r = two*(vel(i+1,j,k,n) - vel(i  ,j,k,n))
                  du_c = -(vel(i-1,j,k,n)+3.d0*vel(i,j,k,n)-4.d0*vel(i+1,j,k,n))/3.d0
                  xslopes(i,j,k,n) = mc_limiter ( du_l, du_c, du_r )

               end if
            end do
         end do
         end do
      end if

      if ( lo(2) == domlo(2) ) then

         j = lo(2)

         do n = 1, 3
         do k = lo(3), hi(3)
            do i = lo(1), hi(1)

               if ( ( bc_jlo_type(i,k,1) == MINF_ ) .or. &
                   ( bc_jlo_type(i,k,1) == NSW_ )  .or. &
                   ( bc_jlo_type(i,k,1) == FSW_ )  .or. &
                   ( bc_jlo_type(i,k,1) == PSW_ )  ) then

                  du_l = two*(vel(i,j  ,k,n) - vel(i,j-1,k,n))
                  du_r = two*(vel(i,j+1,k,n) - vel(i,j  ,k,n))
                  du_c = (vel(i,j+1,k,n)+3.d0*vel(i,j,k,n)-4.d0*vel(i,j-1,k,n))/3.d0
                  yslopes(i,j,k,n) = mc_limiter ( du_l, du_c, du_r )

               end if
            end do
         end do
         end do
      end if

      if ( hi(2) == domhi(2) ) then

         j = hi(2)

         do n = 1, 3
         do k = lo(3), hi(3)
            do i = lo(1), hi(1)

               if ( ( bc_jhi_type(i,k,1) == MINF_ ) .or. &
                   ( bc_jhi_type(i,k,1) == NSW_ )  .or. &
                   ( bc_jhi_type(i,k,1) == FSW_ )  .or. &
                   ( bc_jhi_type(i,k,1) == PSW_ )  ) then

                  du_l = two*(vel(i,j  ,k,n) - vel(i,j-1,k,n))
                  du_r = two*(vel(i,j+1,k,n) - vel(i,j  ,k,n))
                  du_c = -(vel(i,j-1,k,n)+3.d0*vel(i,j,k,n)-4.d0*vel(i,j+1,k,n))/3.d0
                  yslopes(i,j,k,n) = mc_limiter ( du_l, du_c, du_r )

               end if
            end do
         end do
         end do
      end if

      if ( lo(3) == domlo(3) ) then

         k = lo(3)

         do n = 1, 3
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               if ( ( bc_klo_type(i,j,1) == MINF_ ) .or. &
                   ( bc_klo_type(i,j,1) == NSW_ )  .or. &
                   ( bc_klo_type(i,j,1) == FSW_ )  .or. &
                   ( bc_klo_type(i,j,1) == PSW_ )  ) then

                  du_l = two*(vel(i,j,k  ,n) - vel(i,j,k-1,n))
                  du_r = two*(vel(i,j,k+1,n) - vel(i,j,k  ,n))
                  du_c = (vel(i,j,k+1,n)+3.d0*vel(i,j,k,n)-4.d0*vel(i,j,k-1,n))/3.d0
                  zslopes(i,j,k,n) = mc_limiter ( du_l, du_c, du_r )

               end if
            end do
         end do
         end do
      end if

      if ( hi(3) == domhi(3) ) then

         k = hi(3)

         do n = 1, 3
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               if ( ( bc_khi_type(i,j,1) == MINF_ ) .or. &
                   ( bc_khi_type(i,j,1) == NSW_ )  .or. &
                   ( bc_khi_type(i,j,1) == FSW_ )  .or. &
                   ( bc_khi_type(i,j,1) == PSW_ )  ) then

                  du_l = two*(vel(i,j,k  ,n) - vel(i,j,k-1,n))
                  du_r = two*(vel(i,j,k+1,n) - vel(i,j,k  ,n))
                  du_c = -(vel(i,j,k-1,n)+3.d0*vel(i,j,k,n)-4.d0*vel(i,j,k+1,n))/3.d0
                  zslopes(i,j,k,n) = mc_limiter ( du_l, du_c, du_r )

               end if
            end do
         end do
         end do
      end if

   end subroutine compute_slopes

   !
   ! Compute slopes - EB version
   !
   subroutine compute_slopes_eb ( lo, hi, vel, vlo, vhi, &
                                 xslopes, yslopes, zslopes, slo, shi, &
                                 flags, flo, fhi, &
                                 domlo, domhi, &
                                 bc_ilo_type, bc_ihi_type, &
                                 bc_jlo_type, bc_jhi_type, &
                                 bc_klo_type, bc_khi_type, ng ) bind(C)

      use amrex_ebcellflag_module, only: is_covered_cell

      ! Loop bounds
      integer(c_int), intent(in   ) :: lo(3), hi(3)
      integer(c_int), intent(in   ) :: ng

      ! Array bounds
      integer(c_int), intent(in   ) :: vlo(3), vhi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: flo(3), fhi(3)
      ! Grid bounds
      integer(c_int), intent(in   ) :: domlo(3), domhi(3)

      ! BC types
      integer(c_int), intent(in   ) ::  &
           & bc_ilo_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_ihi_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_jlo_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_jhi_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_klo_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
           & bc_khi_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

      ! Arrays
      real(ar),       intent(in   ) ::                      &
           & vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3)

      integer(c_int), intent(in   ) ::                      &
           & flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

      real(ar),       intent(  out) ::                             &
           & xslopes(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3), &
           & yslopes(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3), &
           & zslopes(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      ! Local variables
      integer                       :: i, j, k, n
      real(ar)                      :: du_l, du_c, du_r
      real(ar),    parameter        :: two = 2.0_ar, three2nds = 1.5_ar

      do n = 1, 3
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ! Limiter returns 0 slope if any of the input slopes
                  ! is zero, thus I have just to nullify the correct one sided
                  ! slope in case of neighbor covered cell
                  if ( is_covered_cell(flags(i,j,k)) ) then
                     xslopes(i,j,k,n) = huge(one)
                     yslopes(i,j,k,n) = huge(one)
                     zslopes(i,j,k,n) = huge(one)
                  else
                     ! X direction
                     du_l = two*(vel(i  ,j,k,n) - vel(i-1,j,k,n))
                     du_l = merge( zero, du_l, is_covered_cell(flags(i-1,j,k)) )

                     du_r = two*(vel(i+1,j,k,n) - vel(i  ,j,k,n))
                     du_r = merge( zero, du_r, is_covered_cell(flags(i+1,j,k)) )

                     du_c = half * ( vel(i+1,j,k,n) - vel(i-1,j,k,n) )

                     xslopes(i,j,k,n) = mc_limiter ( du_l, du_c, du_r )

                     ! Y direction
                     du_l = two*(vel(i,j  ,k,n) - vel(i,j-1,k,n))
                     du_l = merge( zero, du_l, is_covered_cell(flags(i,j-1,k)) )

                     du_r = two*(vel(i,j+1,k,n) - vel(i,j  ,k,n))
                     du_r = merge( zero, du_r, is_covered_cell(flags(i,j+1,k)) )

                     du_c = half * ( vel(i,j+1,k,n) - vel(i,j-1,k,n) )

                     yslopes(i,j,k,n) = mc_limiter ( du_l, du_c, du_r )

                     ! z direction
                     du_l = two*(vel(i,j,k  ,n) - vel(i,j,k-1,n))
                     du_l = merge( zero, du_l, is_covered_cell(flags(i,j,k-1)) )

                     du_r = two*(vel(i,j,k+1,n) - vel(i,j,k  ,n))
                     du_r = merge( zero, du_r, is_covered_cell(flags(i,j,k+1)) )

                     du_c = half * ( vel(i,j,k+1,n) - vel(i,j,k-1,n) )

                     zslopes(i,j,k,n) = mc_limiter ( du_l, du_c, du_r )
                  end if
               end do
            end do
         end do
      end do

      !
      ! Compute slopes at boundary where physical BCs are imposed
      !
      if ( lo(1) == domlo(1) ) then

         i = lo(1)

         do n = 1, 3
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)

                  if ( ( bc_ilo_type(j,k,1) == MINF_ ) .or. &
                      ( bc_ilo_type(j,k,1) == NSW_ )  .or. &
                      ( bc_ilo_type(j,k,1) == FSW_ )  .or. &
                      ( bc_ilo_type(j,k,1) == PSW_ )  ) then

                     du_l = two*(vel(i  ,j,k,n) - vel(i-1,j,k,n))
                     du_r = two*(vel(i+1,j,k,n) - vel(i  ,j,k,n))
                     du_c = (vel(i+1,j,k,n)+3.d0*vel(i,j,k,n)-4.d0*vel(i-1,j,k,n))/3.d0
                     xslopes(i,j,k,n) = mc_limiter ( du_l, du_c, du_r )

                  end if
               end do
            end do
         end do
      end if

      if ( hi(1) == domhi(1) ) then

         i = hi(1)

         do n = 1, 3
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)

               if ( ( bc_ihi_type(j,k,1) == MINF_ ) .or. &
                   ( bc_ihi_type(j,k,1) == NSW_ )  .or. &
                   ( bc_ihi_type(j,k,1) == FSW_ )  .or. &
                   ( bc_ihi_type(j,k,1) == PSW_ )  ) then

                  du_l = two*(vel(i  ,j,k,n) - vel(i-1,j,k,n))
                  du_r = two*(vel(i+1,j,k,n) - vel(i  ,j,k,n))
                  du_c = -(vel(i-1,j,k,n)+3.d0*vel(i,j,k,n)-4.d0*vel(i+1,j,k,n))/3.d0
                  xslopes(i,j,k,n) = mc_limiter ( du_l, du_c, du_r )

               end if
            end do
         end do
         end do
      end if

      if ( lo(2) == domlo(2) ) then

         j = lo(2)

         do n = 1, 3
         do k = lo(3), hi(3)
            do i = lo(1), hi(1)

               if ( ( bc_jlo_type(i,k,1) == MINF_ ) .or. &
                   ( bc_jlo_type(i,k,1) == NSW_ )  .or. &
                   ( bc_jlo_type(i,k,1) == FSW_ )  .or. &
                   ( bc_jlo_type(i,k,1) == PSW_ )  ) then

                  du_l = two*(vel(i,j  ,k,n) - vel(i,j-1,k,n))
                  du_r = two*(vel(i,j+1,k,n) - vel(i,j  ,k,n))
                  du_c = (vel(i,j+1,k,n)+3.d0*vel(i,j,k,n)-4.d0*vel(i,j-1,k,n))/3.d0
                  yslopes(i,j,k,n) = mc_limiter ( du_l, du_c, du_r )

               end if
            end do
         end do
         end do
      end if

      if ( hi(2) == domhi(2) ) then

         j = hi(2)

         do n = 1, 3
         do k = lo(3), hi(3)
            do i = lo(1), hi(1)

               if ( ( bc_jhi_type(i,k,1) == MINF_ ) .or. &
                   ( bc_jhi_type(i,k,1) == NSW_ )  .or. &
                   ( bc_jhi_type(i,k,1) == FSW_ )  .or. &
                   ( bc_jhi_type(i,k,1) == PSW_ )  ) then

                  du_l = two*(vel(i,j  ,k,n) - vel(i,j-1,k,n))
                  du_r = two*(vel(i,j+1,k,n) - vel(i,j  ,k,n))
                  du_c = -(vel(i,j-1,k,n)+3.d0*vel(i,j,k,n)-4.d0*vel(i,j+1,k,n))/3.d0
                  yslopes(i,j,k,n) = mc_limiter ( du_l, du_c, du_r )

               end if
            end do
         end do
         end do
      end if

      if ( lo(3) == domlo(3) ) then

         k = lo(3)

         do n = 1, 3
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               if ( ( bc_klo_type(i,j,1) == MINF_ ) .or. &
                   ( bc_klo_type(i,j,1) == NSW_ )  .or. &
                   ( bc_klo_type(i,j,1) == FSW_ )  .or. &
                   ( bc_klo_type(i,j,1) == PSW_ )  ) then

                  du_l = two*(vel(i,j,k  ,n) - vel(i,j,k-1,n))
                  du_r = two*(vel(i,j,k+1,n) - vel(i,j,k  ,n))
                  du_c = (vel(i,j,k+1,n)+3.d0*vel(i,j,k,n)-4.d0*vel(i,j,k-1,n))/3.d0
                  zslopes(i,j,k,n) = mc_limiter ( du_l, du_c, du_r )

               end if
            end do
         end do
         end do
      end if

      if ( hi(3) == domhi(3) ) then

         k = hi(3)

         do n = 1, 3
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               if ( ( bc_khi_type(i,j,1) == MINF_ ) .or. &
                   ( bc_khi_type(i,j,1) == NSW_ )  .or. &
                   ( bc_khi_type(i,j,1) == FSW_ )  .or. &
                   ( bc_khi_type(i,j,1) == PSW_ )  ) then

                  du_l = two*(vel(i,j,k  ,n) - vel(i,j,k-1,n))
                  du_r = two*(vel(i,j,k+1,n) - vel(i,j,k  ,n))
                  du_c = -(vel(i,j,k-1,n)+3.d0*vel(i,j,k,n)-4.d0*vel(i,j,k+1,n))/3.d0
                  zslopes(i,j,k,n) = mc_limiter ( du_l, du_c, du_r )

               end if
            end do
         end do
         end do
      end if

   end subroutine compute_slopes_eb

   !
   ! Monotonized Central (MC) limiter
   !
   function mc_limiter ( dleft, dcenter, dright )  result (slope)

      real(ar), intent(in   ) :: dleft, dcenter, dright
      real(ar)                :: slope
      real(ar), parameter     :: two = 2.0_ar

      slope = min ( abs(dleft), abs(dcenter), abs(dright) )
      slope = merge ( slope, zero, dleft*dright >= zero )
      slope = sign ( one, dcenter ) * slope

   end function mc_limiter

end module slopes_mod
