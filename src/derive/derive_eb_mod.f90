module derive_eb_module

   use amrex_ebcellflag_module, only: is_covered_cell, is_single_valued_cell
   use amrex_fort_module, only: rt => amrex_real
   use iso_c_binding ,    only: c_int

   use constant,          only: zero, half, one, two, my_huge

   implicit none
   private

   public compute_strainrate_eb, compute_vort_eb

contains

   subroutine compute_strainrate_eb(lo, hi,  &
                                    sr, slo, shi,    &
                                    vel, vlo, vhi,    &
                                    flags,    flo,  fhi, &
                                    dx) bind(C)

      ! Loops bounds
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array bounds
      integer(c_int),  intent(in   ) ::  slo(3),  shi(3)
      integer(c_int),  intent(in   ) ::  vlo(3),  vhi(3)
      integer(c_int),  intent(in   ) ::  flo(3),  fhi(3)

      ! Grid
      real(rt), intent(in   ) :: dx(3)

      ! Arrays
      real(rt), intent(in   ) ::   vel(vlo(1):vhi(1), vlo(2):vhi(2), vlo(3):vhi(3), 3)
      real(rt), intent(inout) ::    sr(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3)   )
      integer,  intent(in   ) :: flags(flo(1):fhi(1), flo(2):fhi(2), flo(3):fhi(3)   )

      integer(c_int) :: i, j, k
      real(rt)       :: idx, idy, idz
      real(rt)       :: ux, uy, uz, vx, vy, vz, wx, wy, wz
      real(rt)       :: c0, c1, c2

      ! Coefficients for one-sided difference estimation
      c0 = -1.5d0
      c1 = 2.0d0
      c2 = -0.5d0

      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               if (is_covered_cell(flags(i,j,k))) then

                  sr(i,j,k) = my_huge

               else if (is_single_valued_cell(flags(i,j,k))) then

                  ! Need to check if there are covered cells in neighbours --
                  ! -- if so, use one-sided difference computation (but still quadratic)

                  if (is_covered_cell(flags(i+1,j,k))) then 

                     ! Covered cell to the right, go fish left
                     ux = - (c0 * vel(i,j,k,1) + c1 * vel(i-1,j,k,1) + c2 * vel(i-2,j,k,1)) * idx
                     vx = - (c0 * vel(i,j,k,2) + c1 * vel(i-1,j,k,2) + c2 * vel(i-2,j,k,2)) * idx
                     wx = - (c0 * vel(i,j,k,3) + c1 * vel(i-1,j,k,3) + c2 * vel(i-2,j,k,3)) * idx

                  else if (is_covered_cell(flags(i-1,j,k))) then

                     ! Covered cell to the left, go fish right
                     ux = (c0 * vel(i,j,k,1) + c1 * vel(i+1,j,k,1) + c2 * vel(i+2,j,k,1)) * idx
                     vx = (c0 * vel(i,j,k,2) + c1 * vel(i+1,j,k,2) + c2 * vel(i+2,j,k,2)) * idx
                     wx = (c0 * vel(i,j,k,3) + c1 * vel(i+1,j,k,3) + c2 * vel(i+2,j,k,3)) * idx

                  else

                     ! No covered cells right or left, use standard stencil
                     ux = half * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx
                     vx = half * (vel(i+1,j,k,2) - vel(i-1,j,k,2)) * idx
                     wx = half * (vel(i+1,j,k,3) - vel(i-1,j,k,3)) * idx

                  end if
                                                                   
                  ! Do the same in y-direction 

                  if (is_covered_cell(flags(i,j+1,k))) then

                     uy = - (c0 * vel(i,j,k,1) + c1 * vel(i,j-1,k,1) + c2 * vel(i,j-2,k,1)) * idy
                     vy = - (c0 * vel(i,j,k,2) + c1 * vel(i,j-1,k,2) + c2 * vel(i,j-2,k,2)) * idy
                     wy = - (c0 * vel(i,j,k,3) + c1 * vel(i,j-1,k,3) + c2 * vel(i,j-2,k,3)) * idy

                  else if (is_covered_cell(flags(i,j-1,k))) then

                     uy = (c0 * vel(i,j,k,1) + c1 * vel(i,j+1,k,1) + c2 * vel(i,j+2,k,1)) * idy
                     vy = (c0 * vel(i,j,k,2) + c1 * vel(i,j+1,k,2) + c2 * vel(i,j+2,k,2)) * idy
                     wy = (c0 * vel(i,j,k,3) + c1 * vel(i,j+1,k,3) + c2 * vel(i,j+2,k,3)) * idy
                                                                   
                  else

                     uy = half * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy
                     vy = half * (vel(i,j+1,k,2) - vel(i,j-1,k,2)) * idy
                     wy = half * (vel(i,j+1,k,3) - vel(i,j-1,k,3)) * idy
                                                                   
                  end if

                  ! Do the same in z-direction 

                  if (is_covered_cell(flags(i,j,k+1))) then

                     uz = - (c0 * vel(i,j,k,1) + c1 * vel(i,j,k-1,1) + c2 * vel(i,j,k-2,1)) * idz
                     vz = - (c0 * vel(i,j,k,2) + c1 * vel(i,j,k-1,2) + c2 * vel(i,j,k-2,2)) * idz
                     wz = - (c0 * vel(i,j,k,3) + c1 * vel(i,j,k-1,3) + c2 * vel(i,j,k-2,3)) * idz

                  else if (is_covered_cell(flags(i,j,k-1))) then

                     uz = (c0 * vel(i,j,k,1) + c1 * vel(i,j,k+1,1) + c2 * vel(i,j,k+2,1)) * idz
                     vz = (c0 * vel(i,j,k,2) + c1 * vel(i,j,k+1,2) + c2 * vel(i,j,k+2,2)) * idz
                     wz = (c0 * vel(i,j,k,3) + c1 * vel(i,j,k+1,3) + c2 * vel(i,j,k+2,3)) * idz
                                                                   
                  else

                     uz = half * (vel(i,j,k+1,1) - vel(i,j,k-1,1)) * idz
                     vz = half * (vel(i,j,k+1,2) - vel(i,j,k-1,2)) * idz
                     wz = half * (vel(i,j,k+1,3) - vel(i,j,k-1,3)) * idz
                                                                   
                  end if
                  
                  ! Finally, compute the strain-rate magnitude
                  sr(i,j,k) = sqrt(two * ux**2 + two * vy**2 + two * wz**2 + & 
                                   (uy + vx)**2 + (vz + wy)**2 + (wx + uz)**2)

               else

                  ux = half * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx
                  vx = half * (vel(i+1,j,k,2) - vel(i-1,j,k,2)) * idx
                  wx = half * (vel(i+1,j,k,3) - vel(i-1,j,k,3)) * idx

                  uy = half * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy
                  vy = half * (vel(i,j+1,k,2) - vel(i,j-1,k,2)) * idy
                  wy = half * (vel(i,j+1,k,3) - vel(i,j-1,k,3)) * idy
                                                                   
                  uz = half * (vel(i,j,k+1,1) - vel(i,j,k-1,1)) * idz
                  vz = half * (vel(i,j,k+1,2) - vel(i,j,k-1,2)) * idz
                  wz = half * (vel(i,j,k+1,3) - vel(i,j,k-1,3)) * idz

                  ! Finally, compute the strain-rate magnitude
                  sr(i,j,k) = sqrt(two * ux**2 + two * vy**2 + two * wz**2 + & 
                                   (uy + vx)**2 + (vz + wy)**2 + (wx + uz)**2)

               end if

            end do
         end do
      end do


   end subroutine compute_strainrate_eb

   subroutine compute_vort_eb(lo, hi,  &
                              vort, slo, shi,    &
                              vel, vlo, vhi,    &
                              flags,    flo,  fhi, &
                              dx) bind(C)

      ! Loops bounds
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array bounds
      integer(c_int),  intent(in   ) ::  slo(3),  shi(3)
      integer(c_int),  intent(in   ) ::  vlo(3),  vhi(3)
      integer(c_int),  intent(in   ) ::  flo(3),  fhi(3)

      ! Grid
      real(rt), intent(in   ) :: dx(3)

      ! Arrays
      real(rt), intent(in   ) ::   vel(vlo(1):vhi(1), vlo(2): vhi(2), vlo(3):vhi(3), 3)
      real(rt), intent(inout) ::  vort(slo(1):shi(1), slo(2): shi(2), slo(3):shi(3)   )
      integer,  intent(in   ) :: flags(flo(1):fhi(1), flo(2): fhi(2), flo(3):fhi(3)   )

      integer(c_int) :: i, j, k
      real(rt)       :: idx, idy, idz
      real(rt)       :: uy, uz, vx, vz, wx, wy
      real(rt)       :: c0, c1, c2

      ! Coefficients for one-sided difference estimation
      c0 = -1.5d0
      c1 = 2.0d0
      c2 = -0.5d0

      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               if (is_covered_cell(flags(i,j,k))) then

                  vort(i,j,k) = my_huge

               else if (is_single_valued_cell(flags(i,j,k))) then

                  ! Need to check if there are covered cells in neighbours --
                  ! -- if so, use one-sided difference computation (but still quadratic)

                  if ( is_covered_cell(flags(i+1,j,k)) ) then 

                     ! Covered cell to the right, go fish left
                     vx = - (c0 * vel(i,j,k,2) + c1 * vel(i-1,j,k,2) + c2 * vel(i-2,j,k,2)) * idx
                     wx = - (c0 * vel(i,j,k,3) + c1 * vel(i-1,j,k,3) + c2 * vel(i-2,j,k,3)) * idx

                  else if ( is_covered_cell(flags(i-1,j,k)) ) then 

                     ! Covered cell to the left, go fish right
                     vx = (c0 * vel(i,j,k,2) + c1 * vel(i+1,j,k,2) + c2 * vel(i+2,j,k,2)) * idx
                     wx = (c0 * vel(i,j,k,3) + c1 * vel(i+1,j,k,3) + c2 * vel(i+2,j,k,3)) * idx

                  else

                     ! No covered cells right or left, use standard stencil
                     vx = half * (vel(i+1,j,k,2) - vel(i-1,j,k,2)) * idx
                     wx = half * (vel(i+1,j,k,3) - vel(i-1,j,k,3)) * idx

                  end if
                                                                   
                  ! Do the same in y-direction 

                  if ( is_covered_cell(flags(i,j+1,k)) ) then 

                     uy = - (c0 * vel(i,j,k,1) + c1 * vel(i,j-1,k,1) + c2 * vel(i,j-2,k,1)) * idy
                     wy = - (c0 * vel(i,j,k,3) + c1 * vel(i,j-1,k,3) + c2 * vel(i,j-2,k,3)) * idy

                  else if ( is_covered_cell(flags(i,j-1,k)) ) then 

                     uy = (c0 * vel(i,j,k,1) + c1 * vel(i,j+1,k,1) + c2 * vel(i,j+2,k,1)) * idy
                     wy = (c0 * vel(i,j,k,3) + c1 * vel(i,j+1,k,3) + c2 * vel(i,j+2,k,3)) * idy
                                                                   
                  else

                     uy = half * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy
                     wy = half * (vel(i,j+1,k,3) - vel(i,j-1,k,3)) * idy
                                                                   
                  end if

                  ! Do the same in z-direction 

                  if ( is_covered_cell(flags(i,j,k+1)) ) then 

                     uz = - (c0 * vel(i,j,k,1) + c1 * vel(i,j,k-1,1) + c2 * vel(i,j,k-2,1)) * idz
                     vz = - (c0 * vel(i,j,k,2) + c1 * vel(i,j,k-1,2) + c2 * vel(i,j,k-2,2)) * idz

                  else if ( is_covered_cell(flags(i,j,k-1)) ) then 

                     uz = (c0 * vel(i,j,k,1) + c1 * vel(i,j,k+1,1) + c2 * vel(i,j,k+2,1)) * idz
                     vz = (c0 * vel(i,j,k,2) + c1 * vel(i,j,k+1,2) + c2 * vel(i,j,k+2,2)) * idz
                                                                   
                  else

                     uz = half * (vel(i,j,k+1,1) - vel(i,j,k-1,1)) * idz
                     vz = half * (vel(i,j,k+1,2) - vel(i,j,k-1,2)) * idz
                                                                   
                  end if
                  
                  ! Finally, compute the vorticity magnitude
                  vort(i,j,k) = sqrt((wy-vz)**2 + (uz-wx)**2 + (vx-uy)**2)


               else

                  vx = half * (vel(i+1,j,k,2) - vel(i-1,j,k,2)) * idx
                  wx = half * (vel(i+1,j,k,3) - vel(i-1,j,k,3)) * idx
                  uy = half * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy
                  wy = half * (vel(i,j+1,k,3) - vel(i,j-1,k,3)) * idy
                  uz = (vel(i  ,j  ,k+1,1) - vel(i  ,j  ,k-1,1)) * idz * half
                  vz = (vel(i  ,j  ,k+1,2) - vel(i  ,j  ,k-1,2)) * idz * half
                  
                  vort(i,j,k) = sqrt((wy-vz)**2 + (uz-wx)**2 + (vx-uy)**2)

               end if

            end do
         end do
      end do


   end subroutine compute_vort_eb

   !-----------------------------------------------------------------------!

end module derive_eb_module
