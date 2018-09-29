!
!
!  This module contains the subroutines to perform some of the steps of the
!  MAC projection method.
!
!
module macproj_mod

   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one

   implicit none
   private

   ! Define here the unit vectors
   ! This is used to shift index  based on how the variable is staggered
   ! Check e_x, e_y and e_z in incflo_level.H
   integer(c_int), parameter :: e_i(3,3) = reshape ( [1,0,0,0,1,0,0,0,1], [3,3] )

contains

   subroutine compute_bcoeff_mac ( lo, hi, bcoeff, alo, ahi, &
        u_i, ulo, uhi, ro, slo, shi, dir )  bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: alo(3),ahi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3)

      ! Direction
      integer(c_int), intent(in   ) :: dir

      ! Arrays
      real(ar),       intent(in   ) :: &
           u_i(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           ro(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),       intent(  out) :: &
           bcoeff(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

      integer      :: i, j, k, i0, j0, k0
      real(ar)     :: ro_f

      i0 = e_i(dir,1)
      j0 = e_i(dir,2)
      k0 = e_i(dir,3)

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               ro_f = half * ( ro(i,j,k) + ro(i-i0,j-j0,k-k0) )

               bcoeff(i,j,k) = one / ro_f 

            end do
         end do
      end do

   end subroutine compute_bcoeff_mac
   
   !
   ! Compute the cell-centered divergence of {u,v,w}
   ! 
   subroutine compute_mac_diveu ( lo, hi, diveu, slo, shi, u, ulo, uhi, &
        & v, vlo, vhi, w, wlo, whi, dx )  bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Arrays bounds
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3)
      integer(c_int), intent(in   ) :: vlo(3),vhi(3)
      integer(c_int), intent(in   ) :: wlo(3),whi(3)

      ! Grid 
      real(ar),       intent(in   ) :: dx(3)

      ! Array
      real(ar),       intent(  out) :: &
           diveu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),       intent(in   ) :: &
           u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           v(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           w(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      ! Local variables
      integer  :: i, j, k
      real(ar) :: odx, ody, odz
      real(ar) :: eu_n, eu_s, eu_t, eu_b, eu_e, eu_w

      odx = one / dx(1)
      ody = one / dx(2)
      odz = one / dx(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               eu_e = u(i+1,j,k)
               eu_w = u(i  ,j,k)

               eu_n = v(i,j+1,k)
               eu_s = v(i,j  ,k)

               eu_t = w(i,j,k+1)
               eu_b = w(i,j,k  )

               ! Divergence
               diveu(i,j,k) = (eu_e - eu_w) * odx + (eu_n - eu_s) * ody + &
                    &         (eu_t - eu_b) * odz

            end do
         end do
      end do

   end subroutine compute_mac_diveu


   !
   ! Compute the cell-centered divergence of {u,v,w}
   ! 
   subroutine compute_mac_diveu_eb ( lo, hi, diveu, slo, shi, u, ulo, uhi, &
        & v, vlo, vhi, w, wlo, whi, afracx, axlo, axhi, afracy, aylo, ayhi,    &
        & afracz, azlo, azhi, vfrac, vflo, vfhi, flags, flo, fhi, dx )  bind(C)

      use amrex_ebcellflag_module, only: is_covered_cell
      
      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Arrays bounds
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3)
      integer(c_int), intent(in   ) :: vlo(3),vhi(3)
      integer(c_int), intent(in   ) :: wlo(3),whi(3)
      integer(c_int), intent(in   ) :: axlo(3),axhi(3)
      integer(c_int), intent(in   ) :: aylo(3),ayhi(3)
      integer(c_int), intent(in   ) :: azlo(3),azhi(3)
      integer(c_int), intent(in   ) :: vflo(3),vfhi(3)
      integer(c_int), intent(in   ) :: flo(3), fhi(3)
      
      ! Grid 
      real(ar),       intent(in   ) :: dx(3)

      ! Array
      real(ar),       intent(  out) :: &
           diveu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),       intent(in   ) :: &
           u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           v(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           w(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)), &
           afracx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)), &
           afracy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3)), &
           afracz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3)), &
           vfrac(vflo(1):vfhi(1),vflo(2):vfhi(2),vflo(3):vfhi(3))

      integer(c_int), intent(in   ) ::                      &
           & flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
           
      ! Local variables
      integer  :: i, j, k
      real(ar) :: odx, ody, odz
      real(ar) :: eu_n, eu_s, eu_t, eu_b, eu_e, eu_w

      odx = one / dx(1)
      ody = one / dx(2)
      odz = one / dx(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if ( is_covered_cell(flags(i,j,k)) ) then
                  diveu(i,j,k) = huge(one)
               else                  
                  eu_e = u(i+1,j,k)
                  eu_w = u(i  ,j,k)

                  eu_n = v(i,j+1,k)
                  eu_s = v(i,j  ,k)

                  eu_t = w(i,j,k+1)
                  eu_b = w(i,j,k  )

                  ! Divergence
                  diveu(i,j,k) = ( (eu_e * afracx(i+1,j,k) - eu_w * afracx(i,j,k)) * odx + &
                       &           (eu_n * afracy(i,j+1,k) - eu_s * afracy(i,j,k)) * ody + &
                       &           (eu_t * afracz(i,j,k+1) - eu_b * afracz(i,j,k)) * odz ) &
                       &           / vfrac(i,j,k)
               end if
            end do
         end do
      end do

   end subroutine compute_mac_diveu_eb
   

   !
   ! Computes  u_i = u_i + C * (1/ro) * (dphi/dx_i)
   !
   ! u_i      = i-th component of a staggered vector field u.
   !
   ! ro       = cell centered density field
   !
   ! dphidxi  =  i-th component of staggered gradphi 
   ! 
   ! C        = real constant
   !
   ! dir      = 1,2,3 indicates x, y, z direction respectively.
   !  
   subroutine project_mac_velocity ( lo, hi, u_i, ulo, uhi, &
        & dphidxi, glo, ghi, ro, slo, shi, c, dir ) bind (C)

      ! Loop bounds
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array bounds
      integer(c_int),  intent(in   ) :: ulo(3), uhi(3)
      integer(c_int),  intent(in   ) :: glo(3), ghi(3)
      integer(c_int),  intent(in   ) :: slo(3), shi(3)

      ! Grid and time spacing
      real(ar),        intent(in   ) :: c

      ! Direction
      integer(c_int),  intent(in   ) :: dir

      ! Arrays
      real(ar),        intent(in   ) ::                          &
           &      ro(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
           & dphidxi(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3))

      real(ar),        intent(inout) ::                           &
           & u_i(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      ! Local variables
      integer(c_int)                 :: i, j, k, i0, j0, k0
      real(ar)                       :: oro

      i0 = e_i(dir,1)
      j0 = e_i(dir,2)
      k0 = e_i(dir,3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               oro  = half * ( one/ro(i,j,k) + one/ro(i-i0,j-j0,k-k0) )
               u_i(i,j,k) = u_i(i,j,k) + c * oro * dphidxi(i,j,k)
            end do
         end do
      end do

   end subroutine project_mac_velocity

end module macproj_mod
