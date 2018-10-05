module derive_eb_module

   use amrex_fort_module, only: rt => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one, two
   use amrex_ebcellflag_module, only : get_neighbor_cells_int_single

   implicit none
   private

   real(rt), parameter :: weights(0:2) = [0.d0, 1.d0, 0.5d0]

   public compute_strainrate_eb

contains

   subroutine compute_strainrate_eb(lo, hi,  &
                                    sr, slo, shi,    &
                                    vel, vlo, vhi,    &
                                    flags,    flo,  fhi, &
                                    afrac_x, axlo, axhi, &
                                    afrac_y, aylo, ayhi, &
                                    afrac_z, azlo, azhi, &
                                    cent_x,  cxlo, cxhi, &
                                    cent_y,  cylo, cyhi, &
                                    cent_z,  czlo, czhi, &
                                    vfrac,   vflo, vfhi, &
                                    bcent,    blo,  bhi, &
                                    dx) bind(C)

      ! Loops bounds
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array bounds
      integer(c_int),  intent(in   ) ::  slo(3),  shi(3)
      integer(c_int),  intent(in   ) ::  vlo(3),  vhi(3)
      integer(c_int),  intent(in   ) ::  flo(3),  fhi(3)
      integer(c_int),  intent(in   ) :: axlo(3), axhi(3)
      integer(c_int),  intent(in   ) :: aylo(3), ayhi(3)
      integer(c_int),  intent(in   ) :: azlo(3), azhi(3)
      integer(c_int),  intent(in   ) :: cxlo(3), cxhi(3)
      integer(c_int),  intent(in   ) :: cylo(3), cyhi(3)
      integer(c_int),  intent(in   ) :: czlo(3), czhi(3)
      integer(c_int),  intent(in   ) :: vflo(3), vfhi(3)
      integer(c_int),  intent(in   ) ::  blo(3),  bhi(3)

      ! Grid
      real(rt), intent(in   ) :: dx(3)

      ! Arrays
      real(rt), intent(in   ) ::                                  &
           &     vel(  vlo(1):  vhi(1),  vlo(2):  vhi(2),  vlo(3):  vhi(3),3), &
           & afrac_x( axlo(1): axhi(1), axlo(2): axhi(2), axlo(3): axhi(3)  ), &
           & afrac_y( aylo(1): ayhi(1), aylo(2): ayhi(2), aylo(3): ayhi(3)  ), &
           & afrac_z( azlo(1): azhi(1), azlo(2): azhi(2), azlo(3): azhi(3)  ), &
           &  cent_x( cxlo(1): cxhi(1), cxlo(2): cxhi(2), cxlo(3): cxhi(3),2), &
           &  cent_y( cylo(1): cyhi(1), cylo(2): cyhi(2), cylo(3): cyhi(3),2), &
           &  cent_z( czlo(1): czhi(1), czlo(2): czhi(2), czlo(3): czhi(3),2), &
           &   vfrac( vflo(1): vfhi(1), vflo(2): vfhi(2), vflo(3): vfhi(3)  ), &
           &   bcent(  blo(1):  bhi(1),  blo(2):  bhi(2),  blo(3):  bhi(3),3)

      real(rt),  intent(inout) ::                                 &
         sr(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      integer,  intent(in   ) :: &
         flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

      integer(c_int) :: i, j, k, n
      real(rt)       :: idx, idy, idz
      real(rt)       :: ux, uy, uz, vx, vy, vz, wx, wy, wz

      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

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


   end subroutine compute_strainrate_eb

   !-----------------------------------------------------------------------!

end module derive_eb_module
