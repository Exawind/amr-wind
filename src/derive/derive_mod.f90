module derive_module

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int
   use param, only: half, one, two

   implicit none
   private

   real(rt), parameter :: q4 = one / ( two * two )

contains

   !
   ! Compute the magnitude of the rate-of-strain tensor
   !
   subroutine compute_strainrate(lo, hi, &
                                 sr, slo, shi, &
                                 vel, vlo, vhi, &
                                 dx) bind(C)
      
      integer(c_int), intent(in   ) ::  lo(3), hi(3)
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: vlo(3),vhi(3)

      real(rt),   intent(in   ) :: dx(3)

      real(rt), intent(in   ) :: &
         vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3)

      real(rt),   intent(  out) :: &
         sr(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Local variables
      !-----------------------------------------------
      integer      :: i, j, k
      real(rt) :: idx, idy, idz
      real(rt) :: ux, uy, uz, vx, vy, vz, wx, wy, wz

      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               ux = (vel(i+1,j  ,k  ,1) - vel(i-1,j  ,k  ,1)) * idx
               vx = (vel(i+1,j  ,k  ,2) - vel(i-1,j  ,k  ,2)) * idx
               wx = (vel(i+1,j  ,k  ,3) - vel(i-1,j  ,k  ,3)) * idx
                                                                
               uy = (vel(i  ,j+1,k  ,1) - vel(i  ,j-1,k  ,1)) * idy
               vy = (vel(i  ,j+1,k  ,2) - vel(i  ,j-1,k  ,2)) * idy
               wy = (vel(i  ,j+1,k  ,3) - vel(i  ,j-1,k  ,3)) * idy
                                                                
               uz = (vel(i  ,j  ,k+1,1) - vel(i  ,j  ,k-1,1)) * idz
               vz = (vel(i  ,j  ,k+1,2) - vel(i  ,j  ,k-1,2)) * idz
               wz = (vel(i  ,j  ,k+1,3) - vel(i  ,j  ,k-1,3)) * idz
               
               ! The factor half is included here instead of in each of the above
               sr(i,j,k) = half * &
                  sqrt(two * ux**2 + two * vy**2 + two * wz**2 + &
                       (uy + vx)**2 + (vz + wy)**2 + (wx + uz)**2)

            end do
         end do
      end do

   end subroutine compute_strainrate

   !
   ! Compute the vorticity
   !
   subroutine compute_vort ( lo, hi, vort, slo, shi, vel, vlo, vhi, dx) &
      bind(C, name="compute_vort")

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: vlo(3),vhi(3)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      real(rt),   intent(in   ) :: dx(3)
      real(rt),   intent(  out) :: &
         vort(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(rt), intent(in   ) :: &
         vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3)

      ! Local variables
      !-----------------------------------------------
      integer      :: i, j, k
      real(rt) :: idx, idy, idz
      real(rt) :: uy, uz, vx, vz, wx, wy

      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               uy = idy * ( vel(i,j+1,k,1) - vel(i,j-1,k,1))
               uz = idz * ( vel(i,j,k+1,1) - vel(i,j,k-1,1))

               vx = idz * ( vel(i+1,j,k,2) - vel(i-1,j,k,2))
               vz = idz * ( vel(i,j,k+1,2) - vel(i,j,k-1,2))

               wx = idz * ( vel(i+1,j,k,3) - vel(i-1,j,k,3))
               wy = idy * ( vel(i,j+1,k,3) - vel(i,j-1,k,3))

               ! The factor half is included here instead of in each of the above
               vort(i,j,k) = half * sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)

            end do
         end do
      end do

   end subroutine compute_vort

end module derive_module
