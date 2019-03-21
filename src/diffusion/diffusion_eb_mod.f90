module eb_diffusion_mod

   use amrex_error_module,    only: amrex_abort
   use amrex_fort_module,     only: rt => amrex_real
   use amrex_mempool_module,  only: amrex_allocate, amrex_deallocate
   use iso_c_binding,         only: c_int

   use constant,              only: zero, half, one, two

   implicit none
   private

   real(rt), parameter :: weights(0:2) = [0.d0, 1.d0, 0.5d0]

   public compute_divtau_eb

contains

   subroutine compute_divtau_eb(lo, hi,  &
                                divtau, dlo, dhi,    &
                                vel_in, vinlo, vinhi,&
                                eta, ro, slo, shi,   &
                                flags,    flo,  fhi, &
                                afrac_x, axlo, axhi, &
                                afrac_y, aylo, ayhi, &
                                afrac_z, azlo, azhi, &
                                cent_x,  cxlo, cxhi, &
                                cent_y,  cylo, cyhi, &
                                cent_z,  czlo, czhi, &
                                vfrac,   vflo, vfhi, &
                                bcent,    blo,  bhi, &
                                domlo,  domhi,       &
                                bc_ilo, bc_ihi,      &
                                bc_jlo, bc_jhi,      &
                                bc_klo, bc_khi,      &
                                dx, ng,              &
                                do_explicit_diffusion) bind(C)

      use diffusion_mod, only: fill_vel_diff_bc
      use divop_mod,     only: compute_divop

      ! Loops bounds
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Number of ghost cells
      integer(c_int),  intent(in   ) :: ng

      ! Array bounds
      integer(c_int),  intent(in   ) ::  dlo(3),  dhi(3)
      integer(c_int),  intent(in   ) ::vinlo(3),vinhi(3)
      integer(c_int),  intent(in   ) ::  slo(3),  shi(3)
      integer(c_int),  intent(in   ) ::  flo(3),  fhi(3)
      integer(c_int),  intent(in   ) :: axlo(3), axhi(3)
      integer(c_int),  intent(in   ) :: aylo(3), ayhi(3)
      integer(c_int),  intent(in   ) :: azlo(3), azhi(3)
      integer(c_int),  intent(in   ) :: cxlo(3), cxhi(3)
      integer(c_int),  intent(in   ) :: cylo(3), cyhi(3)
      integer(c_int),  intent(in   ) :: czlo(3), czhi(3)
      integer(c_int),  intent(in   ) :: vflo(3), vfhi(3)
      integer(c_int),  intent(in   ) ::  blo(3),  bhi(3)
      integer(c_int),  intent(in   ) ::domlo(3),domhi(3)

      ! Grid
      real(rt), intent(in   ) :: dx(3)

      ! Arrays
      real(rt), intent(in   ) ::                                  &
           &  vel_in(vinlo(1):vinhi(1),vinlo(2):vinhi(2),vinlo(3):vinhi(3),3), &
           &      ro(  slo(1):  shi(1),  slo(2):  shi(2),  slo(3):  shi(3)  ), &
           &     eta(  slo(1):  shi(1),  slo(2):  shi(2),  slo(3):  shi(3)  ), &
           & afrac_x( axlo(1): axhi(1), axlo(2): axhi(2), axlo(3): axhi(3)  ), &
           & afrac_y( aylo(1): ayhi(1), aylo(2): ayhi(2), aylo(3): ayhi(3)  ), &
           & afrac_z( azlo(1): azhi(1), azlo(2): azhi(2), azlo(3): azhi(3)  ), &
           &  cent_x( cxlo(1): cxhi(1), cxlo(2): cxhi(2), cxlo(3): cxhi(3),2), &
           &  cent_y( cylo(1): cyhi(1), cylo(2): cyhi(2), cylo(3): cyhi(3),2), &
           &  cent_z( czlo(1): czhi(1), czlo(2): czhi(2), czlo(3): czhi(3),2), &
           &   vfrac( vflo(1): vfhi(1), vflo(2): vfhi(2), vflo(3): vfhi(3)  ), &
           &   bcent(  blo(1):  bhi(1),  blo(2):  bhi(2),  blo(3):  bhi(3),3)

      real(rt),  intent(inout) ::                                 &
         divtau(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),3)

      integer,  intent(in   ) :: &
         flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

      ! BC types
      integer, intent(in   ) ::  &
         bc_ilo(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
         bc_ihi(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
         bc_jlo(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
         bc_jhi(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
         bc_klo(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
         bc_khi(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

      ! If true  then we include all the diffusive terms in this explicit result
      ! If false then we include all only the off-diagonal terms here -- we do this
      !     by computing the full tensor then subtracting the diagonal terms
      integer(c_int),  intent(in   ) :: do_explicit_diffusion

      ! Temporary array just to handle bc's
      integer(c_int) :: vlo(3), vhi(3)
      real(rt), dimension(:,:,:,:), pointer, contiguous :: vel

      ! Temporary array to handle viscous fluxes at the cell faces (staggered)
      ! Just reserve space for the tile + 3 ghost layers
      integer, parameter :: nh = 3 ! Number of Halo layers
      real(rt) :: fx(lo(1)-nh:hi(1)+nh+1,lo(2)-nh:hi(2)+nh  ,lo(3)-nh:hi(3)+nh  ,3)
      real(rt) :: fy(lo(1)-nh:hi(1)+nh  ,lo(2)-nh:hi(2)+nh+1,lo(3)-nh:hi(3)+nh  ,3)
      real(rt) :: fz(lo(1)-nh:hi(1)+nh  ,lo(2)-nh:hi(2)+nh  ,lo(3)-nh:hi(3)+nh+1,3)

      integer(c_int) :: i, j, k, n
      real(rt)       :: idx, idy, idz

      ! Check number of ghost cells
      if (ng < 4) call amrex_abort( "compute_divtau_eb(): ng must be >= 4")

      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)

      vlo = lo - ng
      vhi = hi + ng
      call amrex_allocate(vel, vlo(1), vhi(1), vlo(2), vhi(2), vlo(3), vhi(3), 1, 3)

      ! Put values into ghost cells so we can easy take derivatives
      ! TODO: is this necessary? We have filled 
      call fill_vel_diff_bc(vel_in, vinlo, vinhi, vel, lo, hi, domlo, domhi, ng, &
                            bc_ilo, bc_ihi, bc_jlo, bc_jhi, bc_klo, bc_khi)

      ! tau_xx, tau_xy, tau_xz on west faces
      call compute_tau_x(vel, vlo, vhi, eta, slo, shi, &
                         flags, flo, fhi, lo, hi, dx, fx, nh, domlo, domhi, do_explicit_diffusion)

      ! tau_yx, tau_yy, tau_yz on south faces
      call compute_tau_y(vel, vlo, vhi, eta, slo, shi, &
                         flags, flo, fhi, lo, hi, dx, fy, nh, domlo, domhi, do_explicit_diffusion)

      ! tau_zx, tau_zy, tau_zz on bottom faces
      call compute_tau_z(vel, vlo, vhi, eta, slo, shi, &
                         flags, flo, fhi, lo, hi, dx, fz, nh, domlo, domhi, do_explicit_diffusion)

      divop: block
         ! Compute div(tau) with EB algorithm
         integer(c_int)  :: fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3)

         fxlo = lo - nh
         fylo = lo - nh
         fzlo = lo - nh

         fxhi = hi + nh + [1,0,0]
         fyhi = hi + nh + [0,1,0]
         fzhi = hi + nh + [0,0,1]

         call compute_divop( lo, hi, &
                            divtau, dlo, dhi, &
                            vel, vlo, vhi, &
                            fx, fxlo, fxhi, &
                            fy, fylo, fyhi, &
                            fz, fzlo, fzhi, &
                            afrac_x, axlo, axhi, &
                            afrac_y, aylo, ayhi, &
                            afrac_z, azlo, azhi, &
                            cent_x, cxlo, cxhi, &
                            cent_y, cylo, cyhi, &
                            cent_z, czlo, czhi,    &
                            flags, flo, fhi, &
                            vfrac, vflo, vfhi, &
                            bcent, blo, bhi, &
                            domlo, domhi, &
                            dx, ng, eta, do_explicit_diffusion)

      end block divop

      ! Divide by ro
      do n = 1, 3
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  divtau(i,j,k,n) = divtau(i,j,k,n) / ro(i,j,k)
               end do
            end do
         end do
      end do

      call amrex_deallocate(vel)

   end subroutine compute_divtau_eb

   !<$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$>!
   !                             WARNING                                   !
   !<$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$>!
   !                                                                       !
   !  compute_tau_<dir> computes the viscous fluxes at the cell faces.     !
   !  The redistribution algorithm in compute_divop() will ignore all the  !
   !  fluxes outside the domain, EXCEPT when the boundary is periodic.     !
   !  Therefore we always compute the fluxes in the whole ghost region,    !
   !  regardless of the BC. compute_divop() will take care of ignore what  !
   !  is not needed.                                                       !
   !                                                                       !
   !<$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$>!

   !-----------------------------------------------------------------------!
   !-----------------------------------------------------------------------!
   !-----------------------------------------------------------------------!
   !-----------------------------------------------------------------------!
   subroutine compute_tau_x(vel, vlo, vhi, eta, slo, shi, &
                            flag, fglo, fghi, lo, hi, dx, tau_x, ng, domlo, domhi, &
                            do_explicit_diffusion)

      use amrex_ebcellflag_module, only: get_neighbor_cells_int_single

      integer,  intent(in   ) ::  vlo(3),  vhi(3)
      integer,  intent(in   ) ::  slo(3),  shi(3)
      integer,  intent(in   ) :: fglo(3), fghi(3)
      integer,  intent(in   ) ::   lo(3),   hi(3)
      integer,  intent(in   ) ::domlo(3),domhi(3)

      real(rt), intent(in   ) :: dx(3), &
                                 vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3), &
                                 eta(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer,  intent(in   ) :: &
         flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))

      real(rt), intent(  out) :: &
         tau_x(lo(1)-ng:, lo(2)-ng:, lo(3)-ng:,1: )

      integer,  intent(in   ) :: ng ! Number of ghost layer to fill

      integer(c_int),  intent(in   ) :: do_explicit_diffusion

      real(rt) :: dudx, dudy, dudz
      real(rt) :: dvdx, dvdy
      real(rt) :: dwdx,       dwdz

      real(rt) :: wlo, whi
      real(rt) :: idx, idy, idz
      real(rt) :: eta_w

      integer  :: i,j,k
      integer  :: jhip, jhim, jlop, jlom
      integer  :: khip, khim, klop, klom

      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)

      tau_x = zero

      do k = lo(3)-ng, hi(3)+ng
         do j = lo(2)-ng, hi(2)+ng
            do i = lo(1)-ng, hi(1)+ng+1

               dudx = (vel(i,j,k,1) - vel(i-1,j,k,1))*idx
               dvdx = (vel(i,j,k,2) - vel(i-1,j,k,2))*idx
               dwdx = (vel(i,j,k,3) - vel(i-1,j,k,3))*idx

               jhip = j + get_neighbor_cells_int_single(flag(i  ,j,k),0, 1,0)
               jhim = j - get_neighbor_cells_int_single(flag(i  ,j,k),0,-1,0)
               jlop = j + get_neighbor_cells_int_single(flag(i-1,j,k),0, 1,0)
               jlom = j - get_neighbor_cells_int_single(flag(i-1,j,k),0,-1,0)

               if (i == domlo(1)) then
                  jlop = jhip
                  jlom = jhim
               else if (i == domhi(1)+1) then
                  jhip = jlop
                  jhim = jlom
               end if

               whi = weights(jhip-jhim)
               wlo = weights(jlop-jlom)

               dudy = (0.5d0*idy) * &
                  ((vel(i  ,jhip,k,1)-vel(i  ,jhim,k,1))*whi &
                  +(vel(i-1,jlop,k,1)-vel(i-1,jlom,k,1))*wlo)

               dvdy = (0.5d0*idy) * &
                  ((vel(i  ,jhip,k,2)-vel(i  ,jhim,k,2))*whi &
                  +(vel(i-1,jlop,k,2)-vel(i-1,jlom,k,2))*wlo)

               khip = k + get_neighbor_cells_int_single(flag(i  ,j,k),0,0, 1)
               khim = k - get_neighbor_cells_int_single(flag(i  ,j,k),0,0,-1)
               klop = k + get_neighbor_cells_int_single(flag(i-1,j,k),0,0, 1)
               klom = k - get_neighbor_cells_int_single(flag(i-1,j,k),0,0,-1)

               if (i == domlo(1)) then
                  klop = khip
                  klom = khim
               else if (i == domhi(1)+1) then
                  khip = klop
                  khim = klom
               end if

               whi = weights(khip-khim)
               wlo = weights(klop-klom)

               dudz = (0.5d0*idz) * &
                  ((vel(i  ,j,khip,1)-vel(i  ,j,khim,1))*whi &
                  +(vel(i-1,j,klop,1)-vel(i-1,j,klom,1))*wlo)

               dwdz = (0.5d0*idz) * &
                  ((vel(i  ,j,khip,3)-vel(i  ,j,khim,3))*whi &
                  +(vel(i-1,j,klop,3)-vel(i-1,j,klom,3))*wlo)

               eta_w = half * (eta(i,j,k) + eta(i-1,j,k))

               tau_x(i,j,k,1) = eta_w * (dudx + dudx)
               tau_x(i,j,k,2) = eta_w * (dudy + dvdx)
               tau_x(i,j,k,3) = eta_w * (dudz + dwdx)

               if (do_explicit_diffusion .eq. 0) then
                  !
                  ! Subtract diagonal terms of stress tensor, to be obtained through
                  ! implicit solve instead.
                  !
                  tau_x(i,j,k,1) = tau_x(i,j,k,1) - eta_w*dudx
                  tau_x(i,j,k,2) = tau_x(i,j,k,2) - eta_w*dvdx
                  tau_x(i,j,k,3) = tau_x(i,j,k,3) - eta_w*dwdx
               end if

            end do
         end do
      end do

   end subroutine compute_tau_x
   !-----------------------------------------------------------------------!
   !-----------------------------------------------------------------------!
   !-----------------------------------------------------------------------!
   !-----------------------------------------------------------------------!

   subroutine compute_tau_y(vel, vlo, vhi, eta, slo, shi, &
                            flag, fglo, fghi, lo, hi, dx, tau_y, ng, domlo, domhi, &
                            do_explicit_diffusion)

      use amrex_ebcellflag_module, only: get_neighbor_cells_int_single

      integer,  intent(in   ) ::  vlo(3),  vhi(3)
      integer,  intent(in   ) ::  slo(3),  shi(3)
      integer,  intent(in   ) :: fglo(3), fghi(3)
      integer,  intent(in   ) ::   lo(3),   hi(3)
      integer,  intent(in   ) ::domlo(3),domhi(3)

      real(rt), intent(in   ) :: dx(3), &
                                 vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3), &
                                 eta(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer,  intent(in   ) :: &
         flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))

      real(rt), intent(  out) :: &
         tau_y(lo(1)-ng:, lo(2)-ng:, lo(3)-ng:,1: )

      integer,  intent(in   ) :: ng ! Number of ghost layer to fill

      integer(c_int),  intent(in   ) :: do_explicit_diffusion

      real(rt) :: dudx, dudy
      real(rt) :: dvdx, dvdy, dvdz
      real(rt) ::       dwdy, dwdz

      real(rt) :: wlo, whi
      real(rt) :: idx, idy, idz
      real(rt) :: eta_s

      integer  :: i,j,k
      integer  :: ihip, ihim, ilop, ilom
      integer  :: khip, khim, klop, klom

      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)

      tau_y = zero

      do k = lo(3)-ng, hi(3)+ng
         do j = lo(2)-ng, hi(2)+ng+1
            do i = lo(1)-ng, hi(1)+ng

               dudy = (vel(i,j,k,1) - vel(i,j-1,k,1))*idy
               dvdy = (vel(i,j,k,2) - vel(i,j-1,k,2))*idy
               dwdy = (vel(i,j,k,3) - vel(i,j-1,k,3))*idy

               ihip = i + get_neighbor_cells_int_single(flag(i,j  ,k), 1,0,0)
               ihim = i - get_neighbor_cells_int_single(flag(i,j  ,k),-1,0,0)
               ilop = i + get_neighbor_cells_int_single(flag(i,j-1,k), 1,0,0)
               ilom = i - get_neighbor_cells_int_single(flag(i,j-1,k),-1,0,0)

               if (j == domlo(2)) then
                  ilop = ihip
                  ilom = ihim
               else if (j == domhi(2)+1) then
                  ihip = ilop
                  ihim = ilom
               end if

               whi = weights(ihip-ihim)
               wlo = weights(ilop-ilom)

               dudx = (0.5d0*idx) * &
                  ((vel(ihip,j  ,k,1)-vel(ihim,j  ,k,1))*whi &
                  +(vel(ilop,j-1,k,1)-vel(ilom,j-1,k,1))*wlo)

               dvdx = (0.5d0*idx) * &
                  ((vel(ihip,j  ,k,2)-vel(ihim,j  ,k,2))*whi &
                  +(vel(ilop,j-1,k,2)-vel(ilom,j-1,k,2))*wlo)

               khip = k + get_neighbor_cells_int_single(flag(i,j  ,k),0,0, 1)
               khim = k - get_neighbor_cells_int_single(flag(i,j  ,k),0,0,-1)
               klop = k + get_neighbor_cells_int_single(flag(i,j-1,k),0,0, 1)
               klom = k - get_neighbor_cells_int_single(flag(i,j-1,k),0,0,-1)

               if (j == domlo(2)) then
                  klop = khip
                  klom = khim
               else if (j == domhi(2)+1) then
                  khip = klop
                  khim = klom
               end if

               whi = weights(khip-khim)
               wlo = weights(klop-klom)

               dvdz = (0.5d0*idz) * &
                  ((vel(i,j  ,khip,2)-vel(i,j  ,khim,2))*whi &
                  +(vel(i,j-1,klop,2)-vel(i,j-1,klom,2))*wlo)

               dwdz = (0.5d0*idz) * &
                  ((vel(i,j  ,khip,3)-vel(i,j  ,khim,3))*whi &
                  +(vel(i,j-1,klop,3)-vel(i,j-1,klom,3))*wlo)

               eta_s = half * (eta(i,j,k) + eta(i,j-1,k))

               tau_y(i,j,k,1) = eta_s * (dudy + dvdx)
               tau_y(i,j,k,2) = eta_s * (dvdy + dvdy)
               tau_y(i,j,k,3) = eta_s * (dwdy + dvdz)

               if (do_explicit_diffusion .eq. 0) then
                  !
                  ! Subtract diagonal terms of stress tensor, to be obtained through
                  ! implicit solve instead.
                  !
                  tau_y(i,j,k,1) = tau_y(i,j,k,1) - eta_s * dudy
                  tau_y(i,j,k,2) = tau_y(i,j,k,2) - eta_s * dvdy
                  tau_y(i,j,k,3) = tau_y(i,j,k,3) - eta_s * dwdy
               end if

            end do
         end do
      end do

   end subroutine compute_tau_y

   !-----------------------------------------------------------------------!
   !-----------------------------------------------------------------------!
   !-----------------------------------------------------------------------!
   !-----------------------------------------------------------------------!

   subroutine compute_tau_z(vel, vlo, vhi, eta, slo, shi, &
                            flag, fglo, fghi, lo, hi, dx, tau_z, ng, domlo, domhi, & 
                            do_explicit_diffusion)

      use amrex_ebcellflag_module, only: get_neighbor_cells_int_single

      integer,  intent(in   ) ::  vlo(3),  vhi(3)
      integer,  intent(in   ) ::  slo(3),  shi(3)
      integer,  intent(in   ) :: fglo(3), fghi(3)
      integer,  intent(in   ) ::   lo(3),   hi(3)
      integer,  intent(in   ) ::domlo(3),domhi(3)

      real(rt), intent(in   ) :: dx(3), &
                                 vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3), &
                                 eta(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer,  intent(in   ) :: &
         flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))

      real(rt), intent(  out) :: &
         tau_z(lo(1)-ng:, lo(2)-ng:, lo(3)-ng:,1: )

      integer,  intent(in   ) :: ng ! Number of ghost layer to fill

      integer(c_int),  intent(in   ) :: do_explicit_diffusion

      real(rt) :: dudx,       dudz
      real(rt) ::       dvdy, dvdz
      real(rt) :: dwdx, dwdy, dwdz

      real(rt) :: wlo, whi
      real(rt) :: idx, idy, idz
      real(rt) :: eta_b

      integer  :: i,j,k
      integer  :: ihip, ihim, ilop, ilom
      integer  :: jhip, jhim, jlop, jlom

      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)

      tau_z = zero

      do k = lo(3)-ng, hi(3)+ng+1
         do j = lo(2)-ng, hi(2)+ng
            do i = lo(1)-ng, hi(1)+ng

               dudz = (vel(i,j,k,1) - vel(i,j,k-1,1))*idz
               dvdz = (vel(i,j,k,2) - vel(i,j,k-1,2))*idz
               dwdz = (vel(i,j,k,3) - vel(i,j,k-1,3))*idz

               ihip = i + get_neighbor_cells_int_single(flag(i,j,k  ), 1,0,0)
               ihim = i - get_neighbor_cells_int_single(flag(i,j,k  ),-1,0,0)
               ilop = i + get_neighbor_cells_int_single(flag(i,j,k-1), 1,0,0)
               ilom = i - get_neighbor_cells_int_single(flag(i,j,k-1),-1,0,0)

               if (k == domlo(3)) then
                  ilop = ihip
                  ilom = ihim
               else if (k == domhi(3)+1) then
                  ihip = ilop
                  ihim = ilom
               end if

               whi = weights(ihip-ihim)
               wlo = weights(ilop-ilom)

               dudx = (0.5d0*idx) * &
                  ((vel(ihip,j,k  ,1)-vel(ihim,j,k  ,1))*whi &
                  +(vel(ilop,j,k-1,1)-vel(ilom,j,k-1,1))*wlo)

               dwdx = (0.5d0*idx) * &
                  ((vel(ihip,j,k  ,3)-vel(ihim,j,k  ,3))*whi &
                  +(vel(ilop,j,k-1,3)-vel(ilom,j,k-1,3))*wlo)

               jhip = j + get_neighbor_cells_int_single(flag(i,j,k  ),0 ,1,0)
               jhim = j - get_neighbor_cells_int_single(flag(i,j,k  ),0,-1,0)
               jlop = j + get_neighbor_cells_int_single(flag(i,j,k-1),0 ,1,0)
               jlom = j - get_neighbor_cells_int_single(flag(i,j,k-1),0,-1,0)

               if (k == domlo(3)) then
                  jlop = jhip
                  jlom = jhim
               else if (k == domhi(3)+1) then
                  jhip = jlop
                  jhim = jlom
               end if

               whi = weights(jhip-jhim)
               wlo = weights(jlop-jlom)

               dvdy = (0.5d0*idy) * &
                  ((vel(i,jhip,k  ,2)-vel(i,jhim,k  ,2))*whi &
                  +(vel(i,jlop,k-1,2)-vel(i,jlom,k-1,2))*wlo)

               dwdy = (0.5d0*idy) * &
                  ((vel(i,jhip,k  ,3)-vel(i,jhim,k  ,3))*whi &
                  +(vel(i,jlop,k-1,3)-vel(i,jlom,k-1,3))*wlo)

               eta_b = half * (eta(i,j,k) + eta(i,j,k-1))

               tau_z(i,j,k,1) = eta_b * (dudz + dwdx)
               tau_z(i,j,k,2) = eta_b * (dvdz + dwdy)
               tau_z(i,j,k,3) = eta_b * (dwdz + dwdz)

               if (do_explicit_diffusion .eq. 0) then
                  !
                  ! Subtract diagonal terms of stress tensor, to be obtained through
                  ! implicit solve instead.
                  !
                  tau_z(i,j,k,1) = tau_z(i,j,k,1) - eta_b * dudz
                  tau_z(i,j,k,2) = tau_z(i,j,k,2) - eta_b * dvdz
                  tau_z(i,j,k,3) = tau_z(i,j,k,3) - eta_b * dwdz
               end if

            end do
         end do
      end do

   end subroutine compute_tau_z

end module eb_diffusion_mod
