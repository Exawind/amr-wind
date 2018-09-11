module eb_diffusion_mod

   use amrex_fort_module, only: rt => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one, two

   implicit none
   private

   real(rt), parameter :: weights(0:2) = [0.d0, 1.d0, 0.5d0]

   public compute_divtau_eb

 contains

   subroutine compute_divtau_eb ( lo, hi, divtau, dlo, dhi, &
        vel_in, vlo, vhi,          &
        mu, lambda, ro, slo, shi, &
        flag, fglo, fghi,          &
        domlo, domhi,              &
        bc_ilo_type, bc_ihi_type,  &
        bc_jlo_type, bc_jhi_type,  &
        bc_klo_type, bc_khi_type, dx, ng, &
        do_explicit_diffusion) bind(C)

   use diffusion_mod, only: fill_vel_diff_bc

   ! Loops bounds
   integer(c_int),  intent(in   ) :: lo(3),  hi(3)

   ! Number of ghost cells
   integer(c_int),  intent(in   ) :: ng

   ! Array bounds
   integer,  intent(in   ) ::   vlo(3),   vhi(3)
   integer,  intent(in   ) ::   slo(3),   shi(3)
   integer,  intent(in   ) ::   dlo(3),   dhi(3)
   integer,  intent(in   ) :: domlo(3), domhi(3)
   integer,  intent(in   ) ::  fglo(3),  fghi(3)

   ! Grid
   real(rt), intent(in   ) :: dx(3)

   ! Arrays
   real(rt), intent(in   ) ::                                &
        vel_in(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3), &
            ro(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),   &
            mu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)),   &
        lambda(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

   real(rt),  intent(inout) ::                               &
        divtau(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),3)

   integer,  intent(in   ) :: &
        flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))

   ! BC types
   integer, intent(in   ) ::  &
        bc_ilo_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
        bc_ihi_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
        bc_jlo_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
        bc_jhi_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
        bc_klo_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
        bc_khi_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

   ! If true  then we include all the diffusive terms in this explicit result
   ! If false then we include all only the off-diagonal terms here -- we do this
   !     by computing the full tensor then subtracting the diagonal terms
   integer,  intent(in   ) :: do_explicit_diffusion


   ! Temporary array just to handle bc's
   real(rt)   ::  vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3)

   ! Temporaty array to handle div(u) at the nodes
   real(rt) :: tau_x(vlo(1):vhi(1)+1,vlo(2):vhi(2)+1,vlo(3)+1:vhi(3),3)
   real(rt) :: tau_y(vlo(1):vhi(1)+1,vlo(2):vhi(2)+1,vlo(3)+1:vhi(3),3)
   real(rt) :: tau_z(vlo(1):vhi(1)+1,vlo(2):vhi(2)+1,vlo(3)+1:vhi(3),3)

   integer(c_int)  :: i, j, k
   real(rt)        :: idx, idy, idz

   idx = one / dx(1)
   idy = one / dx(2)
   idz = one / dx(3)

   ! Put values into ghost cells so we can easy take derivatives
   call fill_vel_diff_bc(vel_in, vel, vlo, vhi, lo, hi, domlo, domhi, ng, &
                         bc_ilo_type, bc_ihi_type, &
                         bc_jlo_type, bc_jhi_type, &
                         bc_klo_type, bc_khi_type)


   ! tau_xx, tau_xy, tau_xz on west faces
   call compute_tau_x(vel, vlo, vhi, mu, slo, shi, lambda, &
                      flag, fglo, fghi, lo, hi, dx, tau_x)

    ! tau_yx, tau_yy, tau_yz on south faces
   call compute_tau_y(vel, vlo, vhi, mu, slo, shi, lambda, &
                      flag, fglo, fghi, lo, hi, dx, tau_y)

   !  ! tau_zx, tau_zy, tau_zz on bottom faces
   call compute_tau_z(vel, vlo, vhi, mu, slo, shi, lambda, &
                      flag, fglo, fghi, lo, hi, dx, tau_z)


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             divtau(i,j,k,1) =(( tau_x(i+1,j  ,k  ,1) - tau_x(i,j,k,1) ) * idx  + &
                               ( tau_y(i  ,j+1,k  ,1) - tau_y(i,j,k,1) ) * idy  + &
                               ( tau_z(i  ,j  ,k+1,1) - tau_z(i,j,k,1) ) * idz ) / ro(i,j,k)


             divtau(i,j,k,2) =(( tau_x(i+1,j  ,k  ,2) - tau_x(i,j,k,2) ) * idx  + &
                               ( tau_y(i  ,j+1,k  ,2) - tau_y(i,j,k,2) ) * idy  + &
                               ( tau_z(i  ,j  ,k+1,2) - tau_z(i,j,k,2) ) * idz ) / ro(i,j,k)

             divtau(i,j,k,3) =(( tau_x(i+1,j  ,k  ,3) - tau_x(i,j,k,3) ) * idx  + &
                               ( tau_y(i  ,j+1,k  ,3) - tau_y(i,j,k,3) ) * idy  + &
                               ( tau_z(i  ,j  ,k+1,3) - tau_z(i,j,k,3) ) * idz ) / ro(i,j,k)

          end do
       end do
    end do

  end subroutine compute_divtau_eb


!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
   subroutine compute_tau_x(vel, vlo, vhi, mu, slo, shi, lambda, &
                            flag, fglo, fghi, lo, hi, dx, tau_x)

   use amrex_ebcellflag_module, only : get_neighbor_cells_int_single

   integer,  intent(in   ) ::  vlo(3),  vhi(3)
   integer,  intent(in   ) ::  slo(3),  shi(3)
   integer,  intent(in   ) :: fglo(3), fghi(3)
   integer,  intent(in   ) ::   lo(3),   hi(3)

   real(rt), intent(in   ) :: dx(3), &
           vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3), &
            mu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
        lambda(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

   integer,  intent(in   ) :: &
         flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))

   real(rt), intent(  out) :: &
        tau_x(vlo(1):vhi(1)+1,vlo(2):vhi(2)+1,vlo(3)+1:vhi(3),3)

   real(rt) :: dudx, dudy, dudz
   real(rt) :: dvdx, dvdy
   real(rt) :: dwdx,       dwdz

   real(rt) :: wlo, whi
   real(rt) :: idx, idy, idz
   real(rt) :: mu_w, lambda_w

   integer  :: i,j,k
   integer  :: jhip, jhim, jlop, jlom
   integer  :: khip, khim, klop, klom

   idx = one / dx(1)
   idy = one / dx(2)
   idz = one / dx(3)

   do k = lo(3), hi(3)
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)+1

            dudx = (vel(i,j,k,1) - vel(i-1,j,k,1))*idx
            dvdx = (vel(i,j,k,2) - vel(i-1,j,k,2))*idx
            dwdx = (vel(i,j,k,3) - vel(i-1,j,k,3))*idx

            jhip = j + get_neighbor_cells_int_single(flag(i  ,j,k),0, 1,0)
            jhim = j - get_neighbor_cells_int_single(flag(i  ,j,k),0,-1,0)
            jlop = j + get_neighbor_cells_int_single(flag(i-1,j,k),0, 1,0)
            jlom = j - get_neighbor_cells_int_single(flag(i-1,j,k),0,-1,0)

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

            whi = weights(khip-khim)
            wlo = weights(klop-klom)

            dudz = (0.5d0*idz) * &
                   ((vel(i  ,j,khip,1)-vel(i  ,j,khim,1))*whi &
                   +(vel(i-1,j,klop,1)-vel(i-1,j,klom,1))*wlo)

            dwdz = (0.5d0*idz) * &
                   ((vel(i  ,j,khip,3)-vel(i  ,j,khim,3))*whi &
                   +(vel(i-1,j,klop,3)-vel(i-1,j,klom,3))*wlo)

            mu_w     = half * (    mu(i,j,k) +     mu(i-1,j,k))
            lambda_w = half * (lambda(i,j,k) + lambda(i-1,j,k))

            tau_x(i,j,k,1) = mu_w*two*dudx + lambda_w*(dudx + dvdy + dwdz)
            tau_x(i,j,k,2) = mu_w*(dudy + dvdx)
            tau_x(i,j,k,3) = mu_w*(dudz + dwdx)


         end do
      end do
   end do
   end subroutine compute_tau_x
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!

   subroutine compute_tau_y(vel, vlo, vhi, mu, slo, shi, lambda, &
                            flag, fglo, fghi, lo, hi, dx, tau_y)

   use amrex_ebcellflag_module, only : get_neighbor_cells_int_single

   integer,  intent(in   ) ::  vlo(3),  vhi(3)
   integer,  intent(in   ) ::  slo(3),  shi(3)
   integer,  intent(in   ) :: fglo(3), fghi(3)
   integer,  intent(in   ) ::   lo(3),   hi(3)

   real(rt), intent(in   ) :: dx(3), &
           vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3), &
            mu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
        lambda(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

   integer,  intent(in   ) :: &
         flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))

   real(rt), intent(  out) :: &
        tau_y(vlo(1):vhi(1)+1,vlo(2):vhi(2)+1,vlo(3)+1:vhi(3),3)

   real(rt) :: dudx, dudy
   real(rt) :: dvdx, dvdy, dvdz
   real(rt) ::       dwdy, dwdz

   real(rt) :: wlo, whi
   real(rt) :: idx, idy, idz
   real(rt) :: mu_s, lambda_s

   integer  :: i,j,k
   integer  :: ihip, ihim, ilop, ilom
   integer  :: khip, khim, klop, klom

   idx = one / dx(1)
   idy = one / dx(2)
   idz = one / dx(3)

   do k = lo(3), hi(3)
      do j = lo(2), hi(2)+1
         do i = lo(1), hi(1)

            dudy = (vel(i,j,k,1) - vel(i,j-1,k,1))*idy
            dvdy = (vel(i,j,k,2) - vel(i,j-1,k,2))*idy
            dwdy = (vel(i,j,k,3) - vel(i,j-1,k,3))*idy

            ihip = i + get_neighbor_cells_int_single(flag(i,j  ,k), 1,0,0)
            ihim = i - get_neighbor_cells_int_single(flag(i,j  ,k),-1,0,0)
            ilop = i + get_neighbor_cells_int_single(flag(i,j-1,k), 1,0,0)
            ilom = i - get_neighbor_cells_int_single(flag(i,j-1,k),-1,0,0)

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

            whi = weights(khip-khim)
            wlo = weights(klop-klom)

            dvdz = (0.5d0*idz) * &
                 ((vel(i,j  ,khip,2)-vel(i,j  ,khim,2))*whi &
                 +(vel(i,j-1,klop,2)-vel(i,j-1,klom,2))*wlo)

            dwdz = (0.5d0*idz) * &
                 ((vel(i,j  ,khip,3)-vel(i,j  ,khim,3))*whi &
                 +(vel(i,j-1,klop,3)-vel(i,j-1,klom,3))*wlo)

            mu_s     = half * (    mu(i,j,k) +     mu(i,j-1,k))
            lambda_s = half * (lambda(i,j,k) + lambda(i,j-1,k))

            tau_y(i,j,k,1) = mu_s*(dudy + dvdx)
            tau_y(i,j,k,2) = mu_s*(dvdy + dvdy) + lambda_s*(dudx + dvdy + dwdz)
            tau_y(i,j,k,3) = mu_s*(dvdz + dwdy)

         end do
      end do
   end do

 end subroutine compute_tau_y


!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!

   subroutine compute_tau_z(vel, vlo, vhi, mu, slo, shi, lambda, &
                            flag, fglo, fghi, lo, hi, dx, tau_z)

   use amrex_ebcellflag_module, only : get_neighbor_cells_int_single

   integer,  intent(in   ) ::  vlo(3),  vhi(3)
   integer,  intent(in   ) ::  slo(3),  shi(3)
   integer,  intent(in   ) :: fglo(3), fghi(3)
   integer,  intent(in   ) ::   lo(3),   hi(3)

   real(rt), intent(in   ) :: dx(3), &
           vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3), &
            mu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
        lambda(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

   integer,  intent(in   ) :: &
         flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))

   real(rt), intent(  out) :: &
        tau_z(vlo(1):vhi(1)+1,vlo(2):vhi(2)+1,vlo(3)+1:vhi(3),3)

   real(rt) :: dudx,       dudz
   real(rt) ::       dvdy, dvdz
   real(rt) :: dwdx, dwdy, dwdz

   real(rt) :: wlo, whi
   real(rt) :: idx, idy, idz
   real(rt) :: mu_b, lambda_b

   integer  :: i,j,k
   integer  :: ihip, ihim, ilop, ilom
   integer  :: jhip, jhim, jlop, jlom

   idx = one / dx(1)
   idy = one / dx(2)
   idz = one / dx(3)

   do k = lo(3), hi(3)+1
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)

            dudz = (vel(i,j,k,1) - vel(i,j,k-1,1))*idz
            dvdz = (vel(i,j,k,2) - vel(i,j,k-1,2))*idz
            dwdz = (vel(i,j,k,3) - vel(i,j,k-1,3))*idz

            ihip = i + get_neighbor_cells_int_single(flag(i,j,k  ), 1,0,0)
            ihim = i - get_neighbor_cells_int_single(flag(i,j,k  ),-1,0,0)
            ilop = i + get_neighbor_cells_int_single(flag(i,j,k-1), 1,0,0)
            ilom = i - get_neighbor_cells_int_single(flag(i,j,k-1),-1,0,0)

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

            whi = weights(jhip-jhim)
            wlo = weights(jlop-jlom)

            dvdy = (0.5d0*idy) * &
                 ((vel(i,jhip,k  ,2)-vel(i,jhim,k  ,2))*whi &
                 +(vel(i,jlop,k-1,2)-vel(i,jlom,k-1,2))*wlo)

            dwdy = (0.5d0*idy) * &
                 ((vel(i,jhip,k  ,3)-vel(i,jhim,k  ,3))*whi &
                 +(vel(i,jlop,k-1,3)-vel(i,jlom,k-1,3))*wlo)

            mu_b     = half * (    mu(i,j,k) +     mu(i,j,k-1))
            lambda_b = half * (lambda(i,j,k) + lambda(i,j,k-1))

            tau_z(i,j,k,1) = mu_b*(dudz + dwdx)
            tau_z(i,j,k,2) = mu_b*(dvdz + dwdy)
            tau_z(i,j,k,3) = mu_b*(dwdz + dwdz) + lambda_b*(dudx + dvdy + dwdz)

         end do
      end do
   end do

 end subroutine compute_tau_z


  end module eb_diffusion_mod
