module average_cc_to_fc_mod

   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one

   implicit none
   private

contains

   !
   ! Average to faces in chosen direction  -- note we only average the "idir"th
   !    component of cc onto the idir'th face
   !
   subroutine average_cc_to_fc ( lo, hi, fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi,  &
                                cc, slo, shi) bind(C)

      ! Loop bounds (assumed face centered!)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Arrays bounds
      integer(c_int), intent(in   ) :: fxlo(3),fxhi(3)
      integer(c_int), intent(in   ) :: fylo(3),fyhi(3)
      integer(c_int), intent(in   ) :: fzlo(3),fzhi(3)
      integer(c_int), intent(in   ) :: slo(3),shi(3)

      ! Array
      real(ar),       intent(inout) :: &
         fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3)), &
         fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3)), &
         fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3))

      real(ar),       intent(in   ) :: &
         cc(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      ! Local variables
      integer  :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)+1
               fx(i,j,k) = half * ( cc(i-1,j,k,1) + cc(i,j,k,1) )
            end do
         end do
      end do

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)+1
            do i = lo(1), hi(1)
               fy(i,j,k) = half * ( cc(i,j-1,k,2) + cc(i,j,k,2) )
            end do
         end do
      end do

      do k = lo(3), hi(3)+1
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               fz(i,j,k) = half * ( cc(i,j,k-1,3) + cc(i,j,k,3) )
            end do
         end do
      end do

   end subroutine average_cc_to_fc

end module average_cc_to_fc_mod
