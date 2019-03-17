module derive_module

   use amrex_fort_module,  only : rt => amrex_real
   use iso_c_binding,      only: c_int

   use constant,           only: zero, half, one, two

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

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the state
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag        <=  integer tag array
! ::: tag_lo,hi   => index extent of tag array
! ::: state       => state array
! ::: state_lo,hi => index extent of state array
! ::: set         => integer value to tag cell for refinement
! ::: clear       => integer value to untag cell
! ::: lo,hi       => work region we are allowed to change
! ::: dx          => cell size
! ::: problo      => phys loc of lower left corner of prob domain
! ::: time        => problem evolution time
! ::: level       => refinement level of this array
! ::: -----------------------------------------------------------

subroutine state_error(lo, hi, &
                       tag, tlo, thi, &
                       state, slo, shi, &
                       tagval, clearval, &
                       dx, problo, time) bind(C, name="state_error")

  use amrex_fort_module, only : amrex_real
  use iso_c_binding, only : c_char
  implicit none
  
  integer                :: lo(3), hi(3)
  integer                :: slo(3), shi(3)
  integer                :: tlo(3), thi(3)
  character(kind=c_char) :: tag(tlo(1):thi(1), tlo(2):thi(2), tlo(3):thi(3))
  character(kind=c_char) :: tagval, clearval
  real(amrex_real)       :: state(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3))
  real(amrex_real)       :: problo(3), dx(3), time

  integer          :: i, j, k
  real(amrex_real) :: x, y, z

  ! Tag on regions of high phi
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           x = problo(1) + (i + 0.5) * dx(1)
           y = problo(2) + (j + 0.5) * dx(2)
           z = problo(3) + (k + 0.5) * dx(3)
           if ( (abs(z-3.0) < 1.5) .and. (sqrt((x-2.0)**2 + (y-2.0)**2) < 1.0) ) then ! .and. (state(i,j,k) > 0.001) ) then
              print*,"z=",z,", so we are refining"
              tag(i,j,k) = tagval
           endif
        enddo
     enddo
  enddo

end subroutine state_error


end module derive_module
